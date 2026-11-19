# ==============================================================================
# 02. ZASTOSOWANIE KRYTERIÓW DOŚWIADCZALNYCH
#
# Praca magisterska:
# „Planowanie celowanych eksperymentów proteomicznych z wykorzystaniem środowiska R”
# Autor: Daniel Drabik
#
# Skrypt łączy teoretyczne peptydy z etapu 01 z wynikami DIA-NN. Główna analiza
# dotyczy wariantu Trypsyna 1 DPMC po usunięciu peptydów ragged. Kandydaci są
# oceniani kolejno na podstawie obserwowalności, obecności RT w gradiencie
# 30-minutowym, stanu naładowania prekursora, PEP oraz jakości fragmentacji.
#
# Informacje o fragmentach pochodzą z kolumn Fr.0–Fr.11 głównych raportów
# DIA-NN wygenerowanych z opcją --export-quant. Do dalszej analizy przechodzą
# jony y o ładunku 1+ lub 2+, bez utraty neutralnej, z m/z > 300 i Quantity > 0.
# Nie stosuje się osobnego progu dla Fr.Score; parametr ten jest wykorzystywany
# przy porządkowaniu fragmentów i kandydatów.
#
# Na końcu tworzone są listy Top 2 i Top 3 dla całej FASTA oraz panelu 741 białek.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. PAKIETY
# ------------------------------------------------------------------------------

pakiety <- c("tidyverse", "readxl", "openxlsx", "ggplot2", "arrow")

for (p in pakiety) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
  library(p, character.only = TRUE)
}

options(scipen = 999)

# ------------------------------------------------------------------------------
# 1. USTAWIENIA
# ------------------------------------------------------------------------------

# Można podać ścieżki ręcznie. NULL uruchamia automatyczne wyszukiwanie.
folder_etap01 <- NULL
plik_panel <- NULL

# Opcjonalnie można wskazać katalog nadrzędny zawierający raporty DIA-NN.
# Przy NULL skrypt przeszukuje katalog roboczy, katalog skryptu, Downloads i Desktop.
folder_raportow_diann <- NULL

prog_PEP <- 0.01
min_dobrych_fragmentow <- 3L
minimalne_product_mz <- 300
liczba_fragmentow_ocenianych <- 3L
# Zgodnie z zastosowaniem EXPORT-QUANT do projektowania MRM/PRM wybieramy
# 3 fragmenty o najwyższym średnim Fr.Score; nie narzucamy progu Score.
dozwolone_typy_fragmentow <- c("y")
dozwolone_ladunki_fragmentow <- c(1L, 2L)
dozwolone_ladunki_prekursora <- c(2L, 3L)
liczba_slotow_export_quant <- 12L
oczekiwana_liczba_bialek_FASTA <- 20663L
oczekiwana_liczba_bialek_panelu <- 741L

kat_peptydow <- c("0", "1", "2", "3-5", "6-10", "11-20", ">20")
wzorzec_ragged_N <- "^(KK|KR|RK|RR)"
wzorzec_ragged_C <- "(KK|KR|RK|RR)$"
wzorzec_standardu_iRT <- regex("Biognosys[|]iRT-Kit_WR_fusion", ignore_case = TRUE)

# Standardy Biognosys iRT-Kit. Współczynniki kalibracji są wyznaczane ponownie
# bezpośrednio z raportów DIA-NN dla gradientów 30 i 60 min.
standardy_iRT <- tribble(
  ~Peptyd,              ~iRT,
  "LGGNEQVTR",         -24.92,
  "GAGSSEPVTGLDAK",      0.00,
  "VEATFGVDESNAK",      12.39,
  "YILAGVENSK",         19.79,
  "TPVISGGPYEYR",       28.71,
  "DGLDAASYYAPVR",      33.38,
  "ADVTPADFSEWSK",      54.62,
  "GTFIIDPGGVIR",       70.52,
  "GTFIIDPGGVIK",       70.52,
  "TGFIIDPGGVIR",       70.52,
  "LFLQFGAQGSPFLK",    100.00,
  "FLLQFGAQGSPFLK",    100.00
)

# Cztery eksperymenty używane w finalnej analizie. W tej wersji każdy raport
# MUSI zawierać kolumny Fr.0...Fr.11 wygenerowane przez --export-quant.
# Nie korzystamy z *-lib.parquet do oceny jakości fragmentów.
# W razie niejednoznacznych nazw można wpisać pełne ścieżki ręcznie.
pliki_DIANN_reczne <- list(
  "30min_1" = list(report = NULL),
  "30min_2" = list(report = NULL),
  "60min_r1" = list(report = NULL),
  "60min_r2" = list(report = NULL)
)

eksperymenty_DIANN <- tibble(
  Eksperyment = c("30min_1", "30min_2", "60min_r1", "60min_r2"),
  Gradient = c("30min", "30min", "60min", "60min"),
  Report_aliasy = list(
    c("report_2025_30_min_export_quant_added.parquet"),
    c("report_2026_30_min_export_quant_added.parquet"),
    c("report_2025_60_min_r1_export_quant_added.parquet"),
    c("report_2025_60_min_r2_export_quant_added.parquet")
  )
)

# Ustala katalog, z którego uruchomiono skrypt; jest on używany jako punkt odniesienia dla plików wynikowych.
pobierz_katalog_skryptu <- function() {
  ramki <- sys.frames()
  for (i in rev(seq_along(ramki))) {
    ofile <- ramki[[i]]$ofile
    if (!is.null(ofile) && length(ofile) == 1 && nzchar(ofile)) {
      return(dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE)))
    }
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

katalog_skryptu <- pobierz_katalog_skryptu()
korzenie_szukania <- unique(c(
  if (!is.null(folder_raportow_diann)) normalizePath(folder_raportow_diann, winslash = "/", mustWork = FALSE) else NULL,
  normalizePath(getwd(), winslash = "/", mustWork = FALSE),
  katalog_skryptu,
  normalizePath(path.expand("~/Downloads"), winslash = "/", mustWork = FALSE),
  normalizePath(path.expand("~/Desktop"), winslash = "/", mustWork = FALSE)
))
korzenie_szukania <- korzenie_szukania[dir.exists(korzenie_szukania)]

znacznik <- format(Sys.time(), "%Y%m%d_%H%M%S")
folder_wyniki <- file.path(
  katalog_skryptu,
  paste0("WYNIKI_ETAP_02_KRYTERIA_DOSWIADCZALNE_PO_RAGGED_EXPORT_QUANT_HUMAN_ONLY_", znacznik)
)
folder_tabele <- file.path(folder_wyniki, "01_Tabele_Excel")
folder_wykresy <- file.path(folder_wyniki, "02_Wykresy")
folder_kontrola <- file.path(folder_wyniki, "03_Kontrola")
folder_csv <- file.path(folder_wyniki, "04_CSV")
folder_dane_posrednie <- file.path(folder_wyniki, "00_DANE_POSREDNIE_DIANN")
folder_cache <- file.path(katalog_skryptu, "WEJSCIA_ETAP_02_WYPAKOWANE")

for (f in c(folder_wyniki, folder_tabele, folder_wykresy, folder_kontrola, folder_csv, folder_dane_posrednie, folder_cache)) {
  dir.create(f, recursive = TRUE, showWarnings = FALSE)
}

# ------------------------------------------------------------------------------
# 2. WYSZUKIWANIE DANYCH WEJŚCIOWYCH
# ------------------------------------------------------------------------------

# Sprawdza, czy wskazany katalog zawiera komplet wymaganych plików.
ma_pliki <- function(folder, pliki) {
  dir.exists(folder) && all(file.exists(file.path(folder, pliki)))
}

# Jeżeli znaleziono kilka zgodnych wyników, wybiera wersję o najnowszym czasie modyfikacji.
wybierz_najnowszy <- function(sciezki) {
  sciezki <- unique(sciezki[file.exists(sciezki) | dir.exists(sciezki)])
  if (length(sciezki) == 0) return(NULL)
  sciezki[which.max(file.info(sciezki)$mtime)]
}

lista_folderow_do_przeszukania <- function() {
  unique(unlist(lapply(korzenie_szukania, function(root) {
    c(root, list.dirs(root, recursive = TRUE, full.names = TRUE))
  }), use.names = FALSE))
}

lista_zipow_do_przeszukania <- function() {
  unique(unlist(lapply(korzenie_szukania, function(root) {
    list.files(
      root,
      pattern = "[.]zip$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )
  }), use.names = FALSE))
}

znajdz_folder_z_plikami <- function(wymagane, preferowany_wzorzec = NULL) {
  foldery <- lista_folderow_do_przeszukania()
  kandydaci <- foldery[vapply(foldery, ma_pliki, logical(1), pliki = wymagane)]
  if (length(kandydaci) == 0) return(NULL)

  if (!is.null(preferowany_wzorzec)) {
    preferowane <- kandydaci[str_detect(basename(kandydaci), preferowany_wzorzec)]
    if (length(preferowane) > 0) return(wybierz_najnowszy(preferowane))
  }
  wybierz_najnowszy(kandydaci)
}

zip_zawiera_pliki <- function(sciezka_zip, wymagane) {
  listing <- tryCatch(utils::unzip(sciezka_zip, list = TRUE), error = function(e) NULL)
  if (is.null(listing) || nrow(listing) == 0) return(FALSE)
  all(wymagane %in% basename(listing$Name))
}

wypakuj_wymagane_z_zip <- function(sciezka_zip, wymagane, etykieta) {
  listing <- utils::unzip(sciezka_zip, list = TRUE)
  wybrane <- character(0)

  for (plik in wymagane) {
    pasujace <- listing$Name[basename(listing$Name) == plik]
    if (length(pasujace) == 0) stop("W ZIP brakuje pliku: ", plik)
    wybrane <- c(wybrane, pasujace[which.min(nchar(pasujace))])
  }

  folder_docelowy <- file.path(
    folder_cache,
    paste0(etykieta, "_", tools::file_path_sans_ext(basename(sciezka_zip)))
  )
  dir.create(folder_docelowy, recursive = TRUE, showWarnings = FALSE)

  utils::unzip(
    sciezka_zip,
    files = unique(wybrane),
    exdir = folder_docelowy,
    junkpaths = TRUE,
    overwrite = TRUE
  )

  if (!ma_pliki(folder_docelowy, wymagane)) {
    stop("Nie udało się poprawnie wypakować wymaganych plików z: ", sciezka_zip)
  }
  folder_docelowy
}

# Szuka kompletnego zestawu plików w katalogach, a w razie potrzeby także w archiwach ZIP.
znajdz_lub_wypakuj_zestaw <- function(wymagane, preferowany_folder, preferowany_zip, etykieta) {
  folder <- znajdz_folder_z_plikami(wymagane, preferowany_folder)
  if (!is.null(folder)) return(folder)

  zipy <- lista_zipow_do_przeszukania()
  zipy <- zipy[vapply(zipy, zip_zawiera_pliki, logical(1), wymagane = wymagane)]
  if (length(zipy) == 0) return(NULL)

  preferowane <- zipy[str_detect(basename(zipy), preferowany_zip)]
  wybrany_zip <- if (length(preferowane) > 0) wybierz_najnowszy(preferowane) else wybierz_najnowszy(zipy)
  wypakuj_wymagane_z_zip(wybrany_zip, wymagane, etykieta)
}

wymagane_etap01 <- c(
  "E1_peptydy_unikatowe_Trypsyna_0MC.csv",
  "E1_peptydy_unikatowe_Trypsyna_1MC.csv",
  "E1_peptydy_unikatowe_Trypsyna_1MC_PRZED_RAGGED.csv",
  "E1_podsumowanie_4_warianty.csv"
)

if (is.null(folder_etap01)) {
  folder_etap01 <- znajdz_lub_wypakuj_zestaw(
    wymagane_etap01,
    "^WYNIKI_ETAP_00_01_KOREKTY_SKK_HUMAN_ONLY_",
    "^WYNIKI_ETAP_00_01_KOREKTY_SKK_HUMAN_ONLY_.*[.]zip$",
    "ETAP01"
  )
}

# Odnajduje plik z panelem 741 białek.
znajdz_panel <- function() {
  kandydaci <- unique(unlist(lapply(korzenie_szukania, function(root) {
    list.files(
      root,
      pattern = "^ALL_candidates_741.*[.]xlsx$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )
  }), use.names = FALSE))
  wybierz_najnowszy(kandydaci)
}

if (is.null(plik_panel)) plik_panel <- znajdz_panel()

if (is.null(folder_etap01) || !ma_pliki(folder_etap01, wymagane_etap01)) {
  stop("Nie znaleziono finalnego ETAPU 00/01 HUMAN-ONLY. Ustaw folder_etap01 ręcznie.")
}
if (is.null(plik_panel) || !file.exists(plik_panel)) {
  stop("Nie znaleziono pliku ALL_candidates_741*.xlsx. Ustaw plik_panel ręcznie.")
}

folder_etap01 <- normalizePath(folder_etap01, winslash = "/", mustWork = TRUE)
plik_panel <- normalizePath(plik_panel, winslash = "/", mustWork = TRUE)

cat("ETAP 00/01 HUMAN-ONLY:\n", folder_etap01, "\n\n")
cat("Panel 741 białek:\n", plik_panel, "\n\n")
cat("Raporty DIA-NN zostaną odnalezione bezpośrednio i przetworzone w ETAPIE 02.\n\n")

# ------------------------------------------------------------------------------
# 3. FUNKCJE POMOCNICZE
# ------------------------------------------------------------------------------

# Wczytuje pliki CSV zapisane z separatorem średnikowym i przecinkiem dziesiętnym.
czytaj_csv2 <- function(sciezka) {
  readr::read_csv2(
    sciezka,
    show_col_types = FALSE,
    progress = TRUE,
    locale = readr::locale(
      decimal_mark = ",",
      grouping_mark = " ",
      encoding = "UTF-8"
    )
  )
}

# Ujednolica zapis wartości liczbowych pochodzących z plików tekstowych i arkuszy.
na_liczbe <- function(x) {
  if (is.numeric(x)) return(x)
  readr::parse_number(
    as.character(x),
    locale = readr::locale(decimal_mark = ",", grouping_mark = " ")
  )
}

# Konwertuje typowe tekstowe zapisy wartości logicznych na TRUE/FALSE.
na_logic <- function(x) {
  if (is.logical(x)) return(x)
  as.character(x) %in% c("TRUE", "True", "true", "1", "TAK", "Tak", "tak")
}

czy_ragged <- function(peptyd) {
  str_detect(peptyd, wzorzec_ragged_N) | str_detect(peptyd, wzorzec_ragged_C)
}

# Porządkuje kolumny tabel peptydów teoretycznych z etapu 01 przed połączeniem z danymi DIA-NN.
przygotuj_E1 <- function(df) {
  wymagane <- c("Bialko", "Peptyd")
  brak <- setdiff(wymagane, colnames(df))
  if (length(brak) > 0) stop("Tabela E1 nie zawiera kolumn: ", paste(brak, collapse = ", "))
  if (!"Nazwa_bialka" %in% colnames(df)) df$Nazwa_bialka <- NA_character_

  df %>%
    mutate(
      Bialko = as.character(Bialko),
      Peptyd = as.character(Peptyd),
      Nazwa_bialka = as.character(Nazwa_bialka),
      Dlugosc_peptydu = nchar(Peptyd),
      Ragged = czy_ragged(Peptyd)
    ) %>%
    distinct(Bialko, Peptyd, .keep_all = TRUE)
}

# Porządkuje tabelę z obserwacjami DIA-NN i ujednolica typy danych wykorzystywane w dalszych filtrach.
przygotuj_E2 <- function(df) {
  wymagane <- c(
    "Bialko", "Peptyd", "Zmierzony", "RT_sredni", "Najlepszy_PEP",
    "Max_dobrych_fragmentow", "Max_suma_top_fragmentow",
    "Max_srednia_top_fragmentow", "Max_intensywnosc",
    "N_plikow_zmierzony", "Wybrany_ladunek", "Wybrane_mz"
  )
  brak <- setdiff(wymagane, colnames(df))
  if (length(brak) > 0) {
    stop("Tabela E2 nie zawiera wymaganych kolumn: ", paste(brak, collapse = ", "))
  }

  kolumny_numeryczne <- intersect(
    c(
      "Najlepszy_PEP", "Najlepszy_QValue", "RT_sredni", "RT_start_sredni",
      "RT_stop_sredni", "Max_dobrych_fragmentow", "Max_suma_top_fragmentow",
      "Max_srednia_top_fragmentow", "Max_sredni_score_top3", "Max_min_score_top3",
      "Max_udzial_runow_top3", "Max_intensywnosc", "N_plikow_zmierzony",
      "Wybrany_ladunek", "Wybrane_mz", "Wybrane_mz_obserwowane_DIANN",
      "Roznica_mz_DIANN_minus_teoria", "Roznica_mz_ppm"
    ),
    colnames(df)
  )

  df %>%
    mutate(
      Bialko = as.character(Bialko),
      Peptyd = as.character(Peptyd),
      Zmierzony = na_logic(Zmierzony)
    ) %>%
    mutate(across(all_of(kolumny_numeryczne), na_liczbe)) %>%
    mutate(Wybrany_ladunek = suppressWarnings(as.integer(Wybrany_ladunek))) %>%
    distinct(Bialko, Peptyd, .keep_all = TRUE)
}

usun_QValue <- function(df) {
  df %>% select(-matches("QValue|Q[.]Value|Q_Value", ignore.case = TRUE))
}

licz_podsumowanie <- function(df, zbior, etap_id, etap, opis) {
  df <- df %>% distinct(Bialko, Peptyd, .keep_all = TRUE)
  tibble(
    Zbior = zbior,
    Etap_ID = etap_id,
    Etap = etap,
    Opis = opis,
    Liczba_bialek = n_distinct(df$Bialko),
    Liczba_peptydow = nrow(df)
  )
}

# Nadaje kolejność peptydom w obrębie białka i pozostawia maksymalnie n najlepiej ocenionych kandydatów.
ranking_top_n <- function(df, n = 3L) {
  df %>%
    mutate(
      .rank_N_plikow = replace_na(N_plikow_zmierzony, 0),
      .rank_PEP = replace_na(Najlepszy_PEP, 1),
      .rank_fragmenty = replace_na(Max_dobrych_fragmentow, 0),
      .rank_score_fragmentow = replace_na(Max_sredni_score_top3, -Inf),
      .rank_suma_fragmentow = replace_na(Max_suma_top_fragmentow, 0),
      .rank_srednia_fragmentow = replace_na(Max_srednia_top_fragmentow, 0),
      .rank_intensywnosc = replace_na(Max_intensywnosc, 0)
    ) %>%
    group_by(Bialko) %>%
    arrange(
      desc(.rank_N_plikow),
      .rank_PEP,
      desc(.rank_score_fragmentow),
      desc(.rank_suma_fragmentow),
      desc(.rank_srednia_fragmentow),
      desc(.rank_intensywnosc),
      desc(.rank_fragmenty),
      Peptyd,
      .by_group = TRUE
    ) %>%
    mutate(Pozycja_w_bialku = row_number()) %>%
    slice_head(n = n) %>%
    ungroup() %>%
    select(-starts_with(".rank_"))
}

kategoria_liczby_peptydow <- function(n) {
  case_when(
    n == 0 ~ "0",
    n == 1 ~ "1",
    n == 2 ~ "2",
    n >= 3 & n <= 5 ~ "3-5",
    n >= 6 & n <= 10 ~ "6-10",
    n >= 11 & n <= 20 ~ "11-20",
    n > 20 ~ ">20",
    TRUE ~ NA_character_
  )
}

zrob_dystrybucje <- function(df, wszystkie_bialka, zbior, etap_id, etap) {
  wszystkie_bialka <- unique(as.character(wszystkie_bialka))

  zlicz <- df %>%
    distinct(Bialko, Peptyd) %>%
    count(Bialko, name = "Liczba_peptydow")

  rozklad <- tibble(Bialko = wszystkie_bialka) %>%
    left_join(zlicz, by = "Bialko") %>%
    mutate(
      Liczba_peptydow = replace_na(Liczba_peptydow, 0L),
      Kategoria = kategoria_liczby_peptydow(Liczba_peptydow)
    ) %>%
    count(Kategoria, name = "Liczba_bialek") %>%
    mutate(Kategoria = factor(Kategoria, levels = kat_peptydow)) %>%
    complete(
      Kategoria = factor(kat_peptydow, levels = kat_peptydow),
      fill = list(Liczba_bialek = 0L)
    ) %>%
    mutate(
      Zbior = zbior,
      Etap_ID = etap_id,
      Etap = etap,
      Liczba_bialek_razem = length(wszystkie_bialka),
      Procent = round(100 * Liczba_bialek / Liczba_bialek_razem, 2)
    ) %>%
    select(Zbior, Etap_ID, Etap, Kategoria, Liczba_bialek, Procent)

  rozklad
}

zapisz_excel <- function(df, sciezka) {
  openxlsx::write.xlsx(df, file = sciezka, overwrite = TRUE)
}

zapisz_csv2 <- function(df, nazwa) {
  readr::write_excel_csv2(df, file.path(folder_csv, nazwa))
}

theme_prm <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 3),
      plot.subtitle = element_text(size = base_size),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 25, hjust = 1),
      panel.grid.minor = element_blank(),
      legend.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      plot.caption = element_text(size = base_size - 2, hjust = 0)
    )
}

zapisz_wykres <- function(wykres, nazwa, szerokosc = 10, wysokosc = 6) {
  ggsave(file.path(folder_wykresy, paste0(nazwa, ".png")), wykres,
         width = szerokosc, height = wysokosc, dpi = 300, bg = "white")
  ggsave(file.path(folder_wykresy, paste0(nazwa, ".pdf")), wykres,
         width = szerokosc, height = wysokosc, bg = "white")
}

kontrola_globalnej_unikatowosci <- function(df) {
  df %>%
    distinct(Bialko, Peptyd) %>%
    count(Peptyd, name = "Liczba_bialek") %>%
    summarise(Wynik = all(Liczba_bialek == 1L)) %>%
    pull(Wynik)
}

# ------------------------------------------------------------------------------
# 4. WCZYTANIE E1 I PANELU
# ------------------------------------------------------------------------------

# Po tym kroku dostępne są teoretyczne listy peptydów oraz identyfikatory białek panelowych.

sc_E1_0 <- file.path(folder_etap01, "E1_peptydy_unikatowe_Trypsyna_0MC.csv")
sc_E1_1 <- file.path(folder_etap01, "E1_peptydy_unikatowe_Trypsyna_1MC.csv")
sc_E1_1_przed <- file.path(folder_etap01, "E1_peptydy_unikatowe_Trypsyna_1MC_PRZED_RAGGED.csv")
sc_E1_podsum <- file.path(folder_etap01, "E1_podsumowanie_4_warianty.csv")

E1_trypsyna_0DPMC <- czytaj_csv2(sc_E1_0) %>% przygotuj_E1()
E1_trypsyna_1DPMC <- czytaj_csv2(sc_E1_1) %>% przygotuj_E1()
E1_trypsyna_1DPMC_przed <- czytaj_csv2(sc_E1_1_przed) %>% przygotuj_E1()
podsumowanie_E1 <- czytaj_csv2(sc_E1_podsum)

if (!"Bialka_w_FASTA" %in% colnames(podsumowanie_E1)) {
  stop("E1_podsumowanie_4_warianty.csv nie zawiera kolumny Bialka_w_FASTA.")
}
liczba_bialek_FASTA <- as.integer(na_liczbe(podsumowanie_E1$Bialka_w_FASTA[1]))

panel_df <- readxl::read_excel(plik_panel)
kol_id <- case_when(
  "Bialko" %in% colnames(panel_df) ~ "Bialko",
  "ID" %in% colnames(panel_df) ~ "ID",
  "UniProt" %in% colnames(panel_df) ~ "UniProt",
  TRUE ~ NA_character_
)
if (is.na(kol_id)) stop("W pliku panelu nie znaleziono kolumny Bialko, ID ani UniProt.")

# Pobieranie identyfikatorów bez użycia niekwalifikowanej funkcji discard().
# Nazwa discard jest eksportowana także przez pakiet scales, co może maskować
# purrr::discard() i powodować błąd: "comparison (>=) is not possible for language types".
panel_ids <- panel_df[[kol_id]] %>%
  as.character() %>%
  stringr::str_trim()

panel_ids <- panel_ids[!is.na(panel_ids) & nzchar(panel_ids)] %>%
  unique()

# Mianownik FASTA jako pełna lista rekordów Homo sapiens z E1.
# Pełnej listy białek nie ma w tabelach peptydowych, dlatego dla dystrybucji
# wykorzystujemy identyfikatory reprezentowane oraz kontrolowane rekordy zerowe.
# Nazwy sztucznych rekordów zerowych służą wyłącznie do agregacji.
reprezentowane_FASTA <- union(E1_trypsyna_0DPMC$Bialko, E1_trypsyna_1DPMC$Bialko)
brakujace_n <- liczba_bialek_FASTA - length(reprezentowane_FASTA)
if (brakujace_n < 0) stop("Liczba reprezentowanych białek przekracza mianownik FASTA.")
wszystkie_bialka_FASTA_do_dystrybucji <- c(
  reprezentowane_FASTA,
  if (brakujace_n > 0) paste0("__BRAK_REPREZENTACJI_", seq_len(brakujace_n)) else character(0)
)

# ------------------------------------------------------------------------------
# 5. RAPORTY DIA-NN EXPORT-QUANT, FRAGMENTY I GENEROWANIE TABEL E2
# ------------------------------------------------------------------------------

# Każdy raport jest analizowany osobno, a następnie obserwacje są agregowane do poziomu peptyd–prekursor.

# Każdy eksperyment ma jeden główny raport .parquet. Skrypt celowo NIE szuka
# par report + *-lib.parquet, aby nie powtórzyć błędu metodologicznego.
znajdz_report_diann <- function(report_aliasy, etykieta) {
  reczne <- pliki_DIANN_reczne[[etykieta]]

  if (!is.null(reczne) && !is.null(reczne$report)) {
    if (!file.exists(reczne$report)) {
      stop("Nie znaleziono ręcznie wskazanego raportu dla ", etykieta, ": ", reczne$report)
    }
    return(tibble(
      Etykieta = etykieta,
      Folder = normalizePath(dirname(reczne$report), winslash = "/", mustWork = TRUE),
      Plik_Report = normalizePath(reczne$report, winslash = "/", mustWork = TRUE),
      Modyfikacja = file.info(reczne$report)$mtime
    ))
  }

  foldery <- lista_folderow_do_przeszukania()
  rekordy <- lapply(foldery, function(folder) {
    raporty <- file.path(folder, report_aliasy)
    raporty <- raporty[file.exists(raporty)]
    if (length(raporty) == 0) return(NULL)

    indeks <- which.min(match(
      basename(raporty), report_aliasy,
      nomatch = length(report_aliasy) + 1L
    ))
    raport <- raporty[indeks]

    tibble(
      Etykieta = etykieta,
      Folder = normalizePath(folder, winslash = "/", mustWork = TRUE),
      Plik_Report = normalizePath(raport, winslash = "/", mustWork = TRUE),
      Modyfikacja = file.info(raport)$mtime
    )
  })

  kandydaci <- bind_rows(rekordy)
  if (nrow(kandydaci) == 0) {
    stop(
      "Nie znaleziono raportu EXPORT-QUANT dla: ", etykieta,
      ". Oczekiwane nazwy: ", paste(report_aliasy, collapse = " lub "), "."
    )
  }

  kandydaci %>% arrange(desc(Modyfikacja)) %>% slice(1)
}

raporty_DIANN <- pmap_dfr(
  eksperymenty_DIANN,
  function(Eksperyment, Gradient, Report_aliasy) {
    raport <- znajdz_report_diann(Report_aliasy, Eksperyment)
    raport %>% mutate(Eksperyment = Eksperyment, Gradient = Gradient, .before = 1)
  }
)

manifest_DIANN <- raporty_DIANN %>%
  transmute(
    Eksperyment,
    Gradient,
    Folder,
    Plik_Report,
    Rozmiar_report_bajty = file.info(Plik_Report)$size
  )

zapisz_excel(manifest_DIANN, file.path(folder_kontrola, "Kontrola_manifest_plikow_DIANN.xlsx"))
readr::write_excel_csv2(
  manifest_DIANN,
  file.path(folder_dane_posrednie, "00_manifest_uzytych_plikow_DIANN.csv")
)

cat("Użyte raporty DIA-NN EXPORT-QUANT:\n")
print(manifest_DIANN %>% select(Eksperyment, Gradient, Plik_Report))

minimum_bezpieczne <- function(x) {
  x2 <- x[is.finite(x)]
  if (length(x2) == 0) NA_real_ else min(x2)
}

srednia_bezpieczna <- function(x) {
  x2 <- x[is.finite(x)]
  if (length(x2) == 0) NA_real_ else mean(x2)
}

maksimum_bezpieczne <- function(x) {
  x2 <- x[is.finite(x)]
  if (length(x2) == 0) NA_real_ else max(x2)
}

normalizuj_peptyd <- function(x) {
  stringr::str_replace_all(as.character(x), "[^A-Z]", "")
}

# ------------------------------------------------------------------------------
# TEORETYCZNE m/z PREKURSORA DO LIST PRM
# ------------------------------------------------------------------------------
# Masy reszt są zgodne z częścią teoretyczną pipeline'u.
# Dla finalnego m/z prekursora PRM uwzględniana jest stała
# karbamidometylacja cysteiny (UniMod:4), zastosowana w analizie DIA-NN.
masa_reszt_PRM <- c(
  A = 71.037113805,
  R = 156.101111050,
  N = 114.042927470,
  D = 115.026943065,
  C = 103.009184505,
  E = 129.042593135,
  Q = 128.058577540,
  G = 57.021463735,
  H = 137.058911875,
  I = 113.084064015,
  L = 113.084064015,
  K = 128.094963050,
  M = 131.040484645,
  F = 147.068413945,
  P = 97.052763875,
  S = 87.032028435,
  T = 101.047678505,
  W = 186.079312980,
  Y = 163.063328575,
  V = 99.068413945
)

masa_H2O_PRM <- 18.010564684
masa_protonu_PRM <- 1.007276466812
delta_CAM_C <- 57.021463735

oblicz_mase_peptydu_PRM <- function(peptyd) {
  aa <- strsplit(as.character(peptyd), "", fixed = TRUE)[[1]]
  if (length(aa) == 0L || any(!aa %in% names(masa_reszt_PRM))) return(NA_real_)

  masa_bez_modyfikacji <- sum(masa_reszt_PRM[aa]) + masa_H2O_PRM
  liczba_C <- sum(aa == "C")

  masa_bez_modyfikacji + liczba_C * delta_CAM_C
}

oblicz_mz_PRM <- function(peptyd, z) {
  z <- suppressWarnings(as.integer(z))
  if (length(z) != 1L || is.na(z) || z <= 0L) return(NA_real_)

  masa <- oblicz_mase_peptydu_PRM(peptyd)
  if (!is.finite(masa)) return(NA_real_)

  (masa + z * masa_protonu_PRM) / z
}

kolumny_fragmentow_export_quant <- unlist(lapply(
  0:(liczba_slotow_export_quant - 1L),
  function(i) paste0("Fr.", i, c(".Index", ".Id", ".Quantity", ".Score"))
), use.names = FALSE)

# Wczytuje raport DIA-NN, rozpakowuje informacje Fr.* i oblicza metryki fragmentacyjne dla prekursorów.
wczytaj_report_export_quant <- function(plik, eksperyment, gradient) {
  kolumny_podstawowe <- c(
    "Precursor.Id", "Modified.Sequence", "Stripped.Sequence", "Precursor.Charge",
    "Proteotypic", "Protein.Group", "Protein.Ids", "Q.Value", "PEP",
    "Precursor.Quantity", "RT", "RT.Start", "RT.Stop", "Precursor.Mz"
  )
  kolumny <- c(kolumny_podstawowe, kolumny_fragmentow_export_quant)

  df <- arrow::read_parquet(
    plik,
    col_select = dplyr::any_of(kolumny),
    as_data_frame = TRUE
  )

  brak <- setdiff(kolumny, colnames(df))
  if (length(brak) > 0) {
    stop(
      "Raport ", basename(plik),
      " nie zawiera pełnych danych --export-quant. Brak kolumn: ",
      paste(brak, collapse = ", ")
    )
  }

  report <- df %>%
    select(all_of(kolumny_podstawowe)) %>%
    mutate(Eksperyment = eksperyment, Gradient = gradient)

  # Puste Fr.Id oznaczają niewykorzystane sloty spośród maksymalnie 12 pozycji
  # eksportowanych przez DIA-NN. Nie są błędem parsera i nie są kandydatami
  # do dalszej oceny fragmentacyjnej.
  fragmenty_wszystkie <- df %>%
    select(Precursor.Id, all_of(kolumny_fragmentow_export_quant)) %>%
    pivot_longer(
      cols = all_of(kolumny_fragmentow_export_quant),
      names_to = c("Slot", ".value"),
      names_pattern = "^Fr[.]([0-9]+)[.](Index|Id|Quantity|Score)$"
    ) %>%
    transmute(
      Eksperyment = eksperyment,
      Gradient = gradient,
      Precursor.Id = as.character(Precursor.Id),
      Slot = suppressWarnings(as.integer(Slot)),
      Fr.Index = suppressWarnings(as.integer(Index)),
      Fragment.Id = as.character(Id),
      Quantity = suppressWarnings(as.numeric(Quantity)),
      Score = suppressWarnings(as.numeric(Score))
    ) %>%
    mutate(
      Fragment.Type = stringr::str_match(Fragment.Id, "^([A-Za-z]+)")[, 2],
      Fragment.Number = suppressWarnings(as.integer(
        stringr::str_match(Fragment.Id, "^[A-Za-z]+([0-9]+)")[, 2]
      )),
      Fragment.Charge = suppressWarnings(as.integer(
        stringr::str_match(Fragment.Id, "\\^([0-9]+)")[, 2]
      )),
      Fragment.Annotation = stringr::str_match(
        Fragment.Id, "\\^[0-9]+([^/]*)/"
      )[, 2],
      Product.Mz = suppressWarnings(as.numeric(
        stringr::str_match(Fragment.Id, "/([0-9]+(?:[.][0-9]+)?)$")[, 2]
      )),
      Fragment.Annotation = tidyr::replace_na(Fragment.Annotation, ""),
      Pusty_slot = is.na(Fragment.Id) | stringr::str_trim(Fragment.Id) == "",
      Parser_OK = !is.na(Fragment.Type) &
        !is.na(Fragment.Number) &
        !is.na(Fragment.Charge) &
        !is.na(Product.Mz)
    )

  kontrola <- fragmenty_wszystkie %>%
    summarise(
      Eksperyment = eksperyment,
      Gradient = gradient,
      N_wierszy_fragmentow = n(),
      N_prekursorow = n_distinct(Precursor.Id),
      Min_slot = minimum_bezpieczne(Slot),
      Max_slot = maksimum_bezpieczne(Slot),
      N_pustych_slotow_Id = sum(Pusty_slot, na.rm = TRUE),
      N_parser_nieudany = sum(!Pusty_slot & !Parser_OK, na.rm = TRUE),
      N_quantity_NA = sum(is.na(Quantity)),
      N_score_NA = sum(is.na(Score)),
      N_quantity_zero = sum(!is.na(Quantity) & Quantity == 0),
      Score_min = minimum_bezpieczne(Score),
      Score_mediana = stats::median(Score, na.rm = TRUE),
      Score_max = maksimum_bezpieczne(Score)
    )

  # Redukcja pamięci: przed połączeniem czterech runów pozostawiamy tylko
  # fragmenty spełniające kryteria strukturalne, ale zachowujemy Quantity == 0,
  # ponieważ zero ma znaczenie przy obliczaniu średniego Score między runami.
  fragmenty <- fragmenty_wszystkie %>%
    filter(
      !Pusty_slot,
      Parser_OK,
      tolower(Fragment.Type) %in% dozwolone_typy_fragmentow,
      Fragment.Charge %in% dozwolone_ladunki_fragmentow,
      Fragment.Annotation == "",
      Product.Mz > minimalne_product_mz
    ) %>%
    group_by(
      Eksperyment, Gradient, Precursor.Id,
      Fragment.Type, Fragment.Number, Fragment.Charge, Fragment.Annotation
    ) %>%
    summarise(
      Product.Mz = srednia_bezpieczna(Product.Mz),
      Quantity = maksimum_bezpieczne(Quantity),
      Score = srednia_bezpieczna(Score),
      .groups = "drop"
    )

  rm(df, fragmenty_wszystkie)
  gc(verbose = FALSE)

  list(report = report, fragmenty = fragmenty, kontrola = kontrola)
}

cat("Wczytywanie czterech raportów DIA-NN z --export-quant...\n")
wyniki_DIANN <- pmap(
  raporty_DIANN,
  function(Eksperyment, Gradient, Etykieta, Folder, Plik_Report, Modyfikacja) {
    cat("  -> ", Eksperyment, ": ", basename(Plik_Report), "\n", sep = "")
    wczytaj_report_export_quant(Plik_Report, Eksperyment, Gradient)
  }
)

dane_report <- bind_rows(purrr::map(wyniki_DIANN, "report"))
fragmenty_DIANN <- bind_rows(purrr::map(wyniki_DIANN, "fragmenty"))
kontrola_export_quant_run <- bind_rows(purrr::map(wyniki_DIANN, "kontrola"))
rm(wyniki_DIANN)
gc()

# fragmenty_DIANN są już ograniczone do kandydatów strukturalnych w funkcji
# wczytującej, aby nie przechowywać jednocześnie ~12 rekordów na każdy prekursor.
fragmenty_strukturalne <- fragmenty_DIANN

# Liczba pomiarów danego gradientu, w których konkretny prekursor występuje.
prekursor_runy <- dane_report %>%
  distinct(Gradient, Eksperyment, Precursor.Id) %>%
  count(Gradient, Precursor.Id, name = "N_runow_prekursora")

# Agregacja tego samego fragmentu pomiędzy pomiarami referencyjnymi danego
# gradientu. Sumę Score i Quantity dzielimy przez liczbę runów prekursora,
# dlatego brak fragmentu / Quantity == 0 w jednym z runów nie jest ignorowany.
fragmenty_agregat <- fragmenty_strukturalne %>%
  group_by(
    Gradient, Precursor.Id,
    Fragment.Type, Fragment.Number, Fragment.Charge, Fragment.Annotation
  ) %>%
  summarise(
    Product.Mz = srednia_bezpieczna(Product.Mz),
    N_runow_fragmentu = n_distinct(Eksperyment),
    N_runow_z_quantity = sum(!is.na(Quantity) & Quantity > 0),
    Suma_Score = sum(replace_na(Score, 0), na.rm = TRUE),
    Suma_Quantity = sum(replace_na(Quantity, 0), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(prekursor_runy, by = c("Gradient", "Precursor.Id")) %>%
  mutate(
    Sredni_Score = if_else(
      N_runow_prekursora > 0,
      Suma_Score / N_runow_prekursora,
      NA_real_
    ),
    Srednia_Quantity = if_else(
      N_runow_prekursora > 0,
      Suma_Quantity / N_runow_prekursora,
      NA_real_
    ),
    Udzial_runow_z_quantity = if_else(
      N_runow_prekursora > 0,
      N_runow_z_quantity / N_runow_prekursora,
      NA_real_
    ),
    Kwalifikujacy_fragment = N_runow_z_quantity >= 1L
  )

# Liczba wszystkich kwalifikujących się fragmentów prekursora.
liczba_fragmentow_na_prekursor <- fragmenty_agregat %>%
  filter(Kwalifikujacy_fragment) %>%
  count(Gradient, Precursor.Id, name = "Liczba_dobrych_fragmentow")

# Top 3 wg średniego Score. Quantity i powtarzalność służą tylko do
# deterministycznego rozstrzygania remisów.
mapa_top_fragmentow <- fragmenty_agregat %>%
  filter(Kwalifikujacy_fragment) %>%
  group_by(Gradient, Precursor.Id) %>%
  arrange(
    desc(Sredni_Score),
    desc(Udzial_runow_z_quantity),
    desc(Srednia_Quantity),
    Fragment.Type,
    Fragment.Number,
    Fragment.Charge,
    .by_group = TRUE
  ) %>%
  mutate(Pozycja_fragmentu = row_number()) %>%
  filter(Pozycja_fragmentu <= liczba_fragmentow_ocenianych) %>%
  ungroup()

ocena_fragmentow_gradient <- mapa_top_fragmentow %>%
  group_by(Gradient, Precursor.Id) %>%
  summarise(
    N_top_fragmentow = n(),
    Suma_top_fragmentow = sum(Srednia_Quantity, na.rm = TRUE),
    Srednia_top_fragmentow = mean(Srednia_Quantity, na.rm = TRUE),
    Najmocniejszy_fragment = maksimum_bezpieczne(Srednia_Quantity),
    Sredni_score_top3 = mean(Sredni_Score, na.rm = TRUE),
    Minimalny_score_top3 = minimum_bezpieczne(Sredni_Score),
    Sredni_udzial_runow_top3 = mean(Udzial_runow_z_quantity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    liczba_fragmentow_na_prekursor,
    by = c("Gradient", "Precursor.Id")
  )

# Kontrola metody fragmentacyjnej, zapisywana później do pliku kontrolnego.
kontrola_fragmentow_metoda <- tibble(
  Parametr = c(
    "Źródło fragmentów",
    "Typ fragmentu",
    "Ładunki fragmentu",
    "Minimalne m/z",
    "Warunek ilościowy",
    "Próg Fr.Score",
    "Wybór fragmentów do metryk",
    "Minimalna liczba fragmentów"
  ),
  Wartosc = c(
    "Główny report.parquet z --export-quant; Fr.0-Fr.11",
    paste(dozwolone_typy_fragmentow, collapse = ", "),
    paste(dozwolone_ladunki_fragmentow, collapse = ", "),
    paste0("> ", minimalne_product_mz),
    "Quantity > 0 w co najmniej jednym runie danego gradientu",
    "Brak arbitralnego progu; Score służy do rankingu",
    paste0("Top ", liczba_fragmentow_ocenianych, " wg średniego Score między runami"),
    paste0(">= ", min_dobrych_fragmentow)
  )
)

# Dołączamy ocenę fragmentów na poziomie gradient + prekursor do każdego
# rekordu głównego raportu. Nie ma tu danych z *-lib.parquet.
prekursory_zmierzone <- dane_report %>%
  left_join(
    ocena_fragmentow_gradient,
    by = c("Gradient", "Precursor.Id")
  ) %>%
  mutate(
    Peptyd = normalizuj_peptyd(Stripped.Sequence),
    Precursor.Charge = suppressWarnings(as.integer(Precursor.Charge)),
    PEP = suppressWarnings(as.numeric(PEP)),
    Q.Value = suppressWarnings(as.numeric(Q.Value)),
    Precursor.Quantity = suppressWarnings(as.numeric(Precursor.Quantity)),
    RT = suppressWarnings(as.numeric(RT)),
    RT.Start = suppressWarnings(as.numeric(RT.Start)),
    RT.Stop = suppressWarnings(as.numeric(RT.Stop)),
    Precursor.Mz = suppressWarnings(as.numeric(Precursor.Mz))
  ) %>%
  filter(!is.na(Peptyd), Peptyd != "") %>%
  filter(Precursor.Charge %in% dozwolone_ladunki_prekursora)

# Łączy peptydy teoretyczne z obserwacjami DIA-NN i przypisuje parametry wybranego prekursora.
przypisz_rt_z_DIANN <- function(teoretyczne_unikatowe, gradient_docelowy) {
  dane_g <- prekursory_zmierzone %>% filter(Gradient == gradient_docelowy)

  # Jeden rekord na peptyd, stan naładowania i eksperyment.
  na_prekursor_plik <- dane_g %>%
    group_by(Peptyd, Precursor.Charge, Eksperyment) %>%
    summarise(
      RT = srednia_bezpieczna(RT),
      RT_start = srednia_bezpieczna(RT.Start),
      RT_stop = srednia_bezpieczna(RT.Stop),
      PEP = minimum_bezpieczne(PEP),
      QValue = minimum_bezpieczne(Q.Value),
      Intensywnosc = srednia_bezpieczna(Precursor.Quantity),
      Precursor_Mz = srednia_bezpieczna(Precursor.Mz),
      Liczba_dobrych_fragmentow = maksimum_bezpieczne(Liczba_dobrych_fragmentow),
      Suma_top_fragmentow = maksimum_bezpieczne(Suma_top_fragmentow),
      Srednia_top_fragmentow = maksimum_bezpieczne(Srednia_top_fragmentow),
      Najmocniejszy_fragment = maksimum_bezpieczne(Najmocniejszy_fragment),
      Sredni_score_top3 = maksimum_bezpieczne(Sredni_score_top3),
      Minimalny_score_top3 = maksimum_bezpieczne(Minimalny_score_top3),
      Sredni_udzial_runow_top3 = maksimum_bezpieczne(Sredni_udzial_runow_top3),
      .groups = "drop"
    )

  # Agregacja tego samego stanu naładowania pomiędzy runami danego gradientu.
  na_prekursor <- na_prekursor_plik %>%
    group_by(Peptyd, Precursor.Charge) %>%
    summarise(
      N_plikow_zmierzony = n_distinct(Eksperyment),
      RT = srednia_bezpieczna(RT),
      RT_start = srednia_bezpieczna(RT_start),
      RT_stop = srednia_bezpieczna(RT_stop),
      PEP = minimum_bezpieczne(PEP),
      QValue = minimum_bezpieczne(QValue),
      Intensywnosc = maksimum_bezpieczne(Intensywnosc),
      Precursor_Mz = srednia_bezpieczna(Precursor_Mz),
      Liczba_dobrych_fragmentow = maksimum_bezpieczne(Liczba_dobrych_fragmentow),
      Suma_top_fragmentow = maksimum_bezpieczne(Suma_top_fragmentow),
      Srednia_top_fragmentow = maksimum_bezpieczne(Srednia_top_fragmentow),
      Najmocniejszy_fragment = maksimum_bezpieczne(Najmocniejszy_fragment),
      Sredni_score_top3 = maksimum_bezpieczne(Sredni_score_top3),
      Minimalny_score_top3 = maksimum_bezpieczne(Minimalny_score_top3),
      Sredni_udzial_runow_top3 = maksimum_bezpieczne(Sredni_udzial_runow_top3),
      .groups = "drop"
    ) %>%
    mutate(Czy_PEP_OK = !is.na(PEP) & PEP <= prog_PEP) %>%
    group_by(Peptyd) %>%
    mutate(
      N_prekursorow_PEP_OK = sum(Czy_PEP_OK, na.rm = TRUE),
      Czy_ma_prekursor_PEP_OK = N_prekursorow_PEP_OK > 0L
    ) %>%
    ungroup()

  # --------------------------------------------------------------------------
  # KLUCZOWA LOGIKA v3: PEP -> EXPORT-QUANT
  # --------------------------------------------------------------------------
  # Jeżeli dla peptydu istnieje co najmniej jeden stan naładowania z PEP <= 0,01,
  # odrzucamy wszystkie stany niespełniające PEP i dopiero WŚRÓD POZOSTAŁYCH
  # wybieramy najlepszy prekursor na podstawie metryk fragmentacyjnych EXPORT-QUANT.
  # W ten sposób jakość fragmentów nie może zmienić odpowiedzi na pytanie,
  # czy peptyd przechodzi etap PEP.
  najlepszy_po_PEP <- na_prekursor %>%
    filter(Czy_PEP_OK) %>%
    mutate(
      .rank_fragmenty = replace_na(Liczba_dobrych_fragmentow, 0),
      .rank_score = replace_na(Sredni_score_top3, -Inf),
      .rank_suma = replace_na(Suma_top_fragmentow, 0),
      .rank_srednia = replace_na(Srednia_top_fragmentow, 0),
      .rank_najmocniejszy = replace_na(Najmocniejszy_fragment, 0),
      .rank_intensywnosc = replace_na(Intensywnosc, 0),
      .rank_PEP = replace_na(PEP, 1)
    ) %>%
    group_by(Peptyd) %>%
    arrange(
      desc(.rank_fragmenty),
      desc(.rank_score),
      desc(.rank_suma),
      desc(.rank_srednia),
      .rank_PEP,
      desc(.rank_najmocniejszy),
      desc(.rank_intensywnosc),
      desc(N_plikow_zmierzony),
      Precursor.Charge,
      .by_group = TRUE
    ) %>%
    slice_head(n = 1L) %>%
    ungroup() %>%
    select(-starts_with(".rank_"))

  # Dla peptydów, dla których ŻADEN stan naładowania nie spełnia PEP <= 0,01,
  # wybieramy reprezentatywny prekursor wyłącznie po to, aby zachować informacje
  # wcześniejszych etapów (Zmierzony, RT, ładunek). Tutaj najpierw preferujemy
  # najniższy PEP, następnie większą liczbę runów i większą intensywność.
  najlepszy_przed_PEP <- na_prekursor %>%
    filter(!Czy_ma_prekursor_PEP_OK) %>%
    mutate(
      .rank_PEP = replace_na(PEP, 1),
      .rank_intensywnosc = replace_na(Intensywnosc, 0)
    ) %>%
    group_by(Peptyd) %>%
    arrange(
      .rank_PEP,
      desc(N_plikow_zmierzony),
      desc(.rank_intensywnosc),
      Precursor.Charge,
      .by_group = TRUE
    ) %>%
    slice_head(n = 1L) %>%
    ungroup() %>%
    select(-starts_with(".rank_"))

  najlepszy_prekursor <- bind_rows(
    najlepszy_po_PEP,
    najlepszy_przed_PEP
  ) %>%
    distinct(Peptyd, .keep_all = TRUE) %>%
    rename(
      RT_sredni = RT,
      RT_start_sredni = RT_start,
      RT_stop_sredni = RT_stop,
      Najlepszy_PEP = PEP,
      Najlepszy_QValue = QValue,
      Max_dobrych_fragmentow = Liczba_dobrych_fragmentow,
      Max_suma_top_fragmentow = Suma_top_fragmentow,
      Max_srednia_top_fragmentow = Srednia_top_fragmentow,
      Max_sredni_score_top3 = Sredni_score_top3,
      Max_min_score_top3 = Minimalny_score_top3,
      Max_udzial_runow_top3 = Sredni_udzial_runow_top3,
      Max_intensywnosc = Intensywnosc,
      Wybrany_ladunek = Precursor.Charge,
      Wybrane_mz_obserwowane_DIANN = Precursor_Mz
    ) %>%
    rowwise() %>%
    mutate(
      Wybrane_mz = round(oblicz_mz_PRM(Peptyd, Wybrany_ladunek), 6),
      Roznica_mz_DIANN_minus_teoria = if_else(
        is.finite(Wybrane_mz_obserwowane_DIANN) & is.finite(Wybrane_mz),
        Wybrane_mz_obserwowane_DIANN - Wybrane_mz,
        NA_real_
      ),
      Roznica_mz_ppm = if_else(
        is.finite(Wybrane_mz_obserwowane_DIANN) & is.finite(Wybrane_mz) & Wybrane_mz > 0,
        1e6 * (Wybrane_mz_obserwowane_DIANN - Wybrane_mz) / Wybrane_mz,
        NA_real_
      ),
      Zrodlo_mz = "teoretyczne_z_sekwencji_i_ladunku_CAM_C"
    ) %>%
    ungroup()

  teoretyczne_unikatowe %>%
    left_join(najlepszy_prekursor, by = "Peptyd") %>%
    mutate(
      Zmierzony = !is.na(RT_sredni),
      Czy_PEP_ponizej_0_01 = !is.na(Najlepszy_PEP) & Najlepszy_PEP <= prog_PEP
    )
}

# Tabele E2 są odtwarzane od początku dla aktualnego E1 HUMAN-ONLY.
cat("Generowanie E2 Trypsyna 0 DPMC, 30 min...\n")
E2_0_30_raw <- przypisz_rt_z_DIANN(E1_trypsyna_0DPMC, "30min") %>% przygotuj_E2()
cat("Generowanie E2 Trypsyna 0 DPMC, 60 min...\n")
E2_0_60_raw <- przypisz_rt_z_DIANN(E1_trypsyna_0DPMC, "60min") %>% przygotuj_E2()

# Wariant 1 DPMC jest najpierw generowany przed filtrem ragged, a następnie
# ograniczany do aktualnego E1 po filtrze. Zapewnia to jedną wspólną ocenę DIA-NN
# i unika dwukrotnego przetwarzania tego samego dużego zbioru.
cat("Generowanie E2 Trypsyna 1 DPMC przed ragged, 30 min...\n")
E2_1_30_przed_ragged <- przypisz_rt_z_DIANN(E1_trypsyna_1DPMC_przed, "30min") %>% przygotuj_E2()
cat("Generowanie E2 Trypsyna 1 DPMC przed ragged, 60 min...\n")
E2_1_60_przed_ragged <- przypisz_rt_z_DIANN(E1_trypsyna_1DPMC_przed, "60min") %>% przygotuj_E2()

E2_1_30_raw <- E2_1_30_przed_ragged %>%
  semi_join(E1_trypsyna_1DPMC %>% distinct(Bialko, Peptyd), by = c("Bialko", "Peptyd"))
E2_1_60_raw <- E2_1_60_przed_ragged %>%
  semi_join(E1_trypsyna_1DPMC %>% distinct(Bialko, Peptyd), by = c("Bialko", "Peptyd"))

# Kalibracja iRT jest generowana w tym samym etapie z tych samych raportów DIA-NN.
rt_standardy <- function(gradient_docelowy) {
  prekursory_zmierzone %>%
    filter(Gradient == gradient_docelowy, Peptyd %in% standardy_iRT$Peptyd) %>%
    group_by(Peptyd) %>%
    summarise(RT = srednia_bezpieczna(RT), .groups = "drop") %>%
    inner_join(standardy_iRT, by = "Peptyd") %>%
    arrange(iRT, Peptyd)
}

kalib_30 <- rt_standardy("30min")
kalib_60 <- rt_standardy("60min")

if (nrow(kalib_30) < 2L || nrow(kalib_60) < 2L) {
  stop("Za mało standardów iRT do kalibracji. Sprawdź raporty DIA-NN 30 i 60 min.")
}

model_30 <- lm(iRT ~ RT, data = kalib_30)
model_60 <- lm(iRT ~ RT, data = kalib_60)

kalibracja_iRT <- tibble(
  Gradient = c("30min", "60min"),
  Nachylenie_a = c(coef(model_30)[["RT"]], coef(model_60)[["RT"]]),
  Przeciecie_b = c(coef(model_30)[["(Intercept)"]], coef(model_60)[["(Intercept)"]]),
  R2 = c(summary(model_30)$r.squared, summary(model_60)$r.squared),
  N_standardow = c(nrow(kalib_30), nrow(kalib_60))
)

punkty_kalibracyjne <- bind_rows(
  kalib_30 %>% mutate(Gradient = "30min"),
  kalib_60 %>% mutate(Gradient = "60min")
) %>% select(Gradient, Peptyd, iRT, RT)

# Eksport danych pośrednich używanych przez ETAP 05.
readr::write_excel_csv2(E2_0_30_raw, file.path(folder_dane_posrednie, "E2_Trypsyna_0MC_RT_30min.csv"))
readr::write_excel_csv2(E2_0_60_raw, file.path(folder_dane_posrednie, "E2_Trypsyna_0MC_RT_60min.csv"))
readr::write_excel_csv2(E2_1_30_raw, file.path(folder_dane_posrednie, "E2_Trypsyna_1MC_RT_30min.csv"))
readr::write_excel_csv2(E2_1_60_raw, file.path(folder_dane_posrednie, "E2_Trypsyna_1MC_RT_60min.csv"))
readr::write_excel_csv2(kalibracja_iRT, file.path(folder_dane_posrednie, "E4_kalibracja_iRT.csv"))
readr::write_excel_csv2(punkty_kalibracyjne, file.path(folder_dane_posrednie, "E4_punkty_kalibracyjne_iRT.csv"))

klucze_E1_0 <- E1_trypsyna_0DPMC %>% distinct(Bialko, Peptyd)
klucze_E1_1 <- E1_trypsyna_1DPMC %>% distinct(Bialko, Peptyd)
klucze_E1_1_przed <- E1_trypsyna_1DPMC_przed %>% distinct(Bialko, Peptyd)
klucze_ragged_usuniete <- anti_join(klucze_E1_1_przed, klucze_E1_1, by = c("Bialko", "Peptyd"))

kontrola_pokrycia_E2 <- tibble()

sprawdz_i_ogranicz_E2 <- function(E2_raw, klucze, etykieta) {
  braki <- anti_join(klucze, E2_raw, by = c("Bialko", "Peptyd"))
  ograniczone <- E2_raw %>% semi_join(klucze, by = c("Bialko", "Peptyd"))

  wiersz <- tibble(
    Tabela = etykieta,
    Pary_E1 = nrow(klucze),
    Pary_E2_po_ograniczeniu = nrow(ograniczone),
    Brakujace_pary_E1_w_E2 = nrow(braki),
    Pokrycie_proc = round(100 * nrow(ograniczone) / nrow(klucze), 4)
  )
  kontrola_pokrycia_E2 <<- bind_rows(kontrola_pokrycia_E2, wiersz)

  if (nrow(braki) > 0) {
    zapisz_excel(
      braki,
      file.path(folder_kontrola, paste0("BLAD_brakujace_pary_", etykieta, ".xlsx"))
    )
    stop("Tabela ", etykieta, " nie pokrywa całego aktualnego E1.")
  }
  ograniczone
}

E2_0_30 <- sprawdz_i_ogranicz_E2(E2_0_30_raw, klucze_E1_0, "Trypsyna_0DPMC_30min")
E2_0_60 <- sprawdz_i_ogranicz_E2(E2_0_60_raw, klucze_E1_0, "Trypsyna_0DPMC_60min")
E2_1_30 <- sprawdz_i_ogranicz_E2(E2_1_30_raw, klucze_E1_1, "Trypsyna_1DPMC_30min_PO_RAGGED")
E2_1_60 <- sprawdz_i_ogranicz_E2(E2_1_60_raw, klucze_E1_1, "Trypsyna_1DPMC_60min_PO_RAGGED")

# Aliasy używane w dalszej części skryptu dla głównego wariantu 1 DPMC.
E2_30 <- E2_1_30
E2_60 <- E2_1_60

# Rekordy usunięte przez ragged są potrzebne wyłącznie do porównania przed/po.
E2_30_ragged_usuniete <- E2_1_30_przed_ragged %>%
  semi_join(klucze_ragged_usuniete, by = c("Bialko", "Peptyd"))
E2_60_ragged_usuniete <- E2_1_60_przed_ragged %>%
  semi_join(klucze_ragged_usuniete, by = c("Bialko", "Peptyd"))

# Kontrola informacyjna: ile rekordów usuniętych jako ragged było wykrywanych
# w danych DIA-NN przed zastosowaniem filtra.
kontrola_ragged_E2 <- bind_rows(
  E2_1_30_przed_ragged %>%
    semi_join(klucze_ragged_usuniete, by = c("Bialko", "Peptyd")) %>%
    summarise(Gradient = "30min", Ragged_usuniete = n(), Ragged_zmierzone = sum(Zmierzony, na.rm = TRUE)),
  E2_1_60_przed_ragged %>%
    semi_join(klucze_ragged_usuniete, by = c("Bialko", "Peptyd")) %>%
    summarise(Gradient = "60min", Ragged_usuniete = n(), Ragged_zmierzone = sum(Zmierzony, na.rm = TRUE))
)

rm(dane_report, fragmenty_DIANN, fragmenty_strukturalne, fragmenty_agregat, liczba_fragmentow_na_prekursor)
gc()

# ------------------------------------------------------------------------------
# 6. PORÓWNANIE 0 I 1 DPMC DLA GRADIENTÓW 30 I 60 MIN
# ------------------------------------------------------------------------------

# Tworzy kolejne poziomy selekcji dla jednego gradientu, aby policzyć wpływ poszczególnych kryteriów.
zbuduj_zbiory_gradient <- function(start_df, E2_df, panel_ids) {
  start_all <- start_df %>% distinct(Bialko, Peptyd, .keep_all = TRUE)
  start_panel <- start_all %>% filter(Bialko %in% panel_ids)

  measured_all <- E2_df %>%
    filter(Zmierzony) %>%
    distinct(Bialko, Peptyd, .keep_all = TRUE)
  measured_panel <- measured_all %>% filter(Bialko %in% panel_ids)

  rt_all <- measured_all %>%
    filter(!is.na(RT_sredni), Wybrany_ladunek %in% dozwolone_ladunki_prekursora)
  rt_panel <- rt_all %>% filter(Bialko %in% panel_ids)

  pep_all <- rt_all %>%
    filter(!is.na(Najlepszy_PEP), Najlepszy_PEP <= prog_PEP)
  pep_panel <- pep_all %>% filter(Bialko %in% panel_ids)

  frag_all <- pep_all %>%
    filter(!is.na(Max_dobrych_fragmentow), Max_dobrych_fragmentow >= min_dobrych_fragmentow)
  frag_panel <- frag_all %>% filter(Bialko %in% panel_ids)

  list(
    start_all = start_all, start_panel = start_panel,
    measured_all = measured_all, measured_panel = measured_panel,
    rt_all = rt_all, rt_panel = rt_panel,
    pep_all = pep_all, pep_panel = pep_panel,
    frag_all = frag_all, frag_panel = frag_panel
  )
}

podsumuj_wariant_gradient <- function(z, wariant, wariant_id, gradient) {
  zrob <- function(zbior, start, measured, rt, pep, frag) {
    bind_rows(
      licz_podsumowanie(start, zbior, 1L, "Peptydy teoretyczne", "Peptydy unikatowe z aktualnego E1 HUMAN-ONLY."),
      licz_podsumowanie(measured, zbior, 2L, "Zmierzone w danym gradiencie", "Peptydy oznaczone jako zmierzone w tabeli E2."),
      licz_podsumowanie(rt, zbior, 3L, "Obserwowany RT i ładunek 2+/3+", "Peptydy z obserwowanym RT i wybranym prekursorem 2+ lub 3+."),
      licz_podsumowanie(pep, zbior, 4L, "PEP <= 0,01", "Peptydy spełniające próg PEP."),
      licz_podsumowanie(frag, zbior, 5L, "PEP <= 0,01 i >=3 fragmenty", "Peptydy spełniające PEP i kryterium fragmentacyjne.")
    ) %>%
      mutate(
        Wariant = wariant,
        Wariant_ID = wariant_id,
        Gradient = gradient,
        Procent_bialek_wzgledem_startu = round(100 * Liczba_bialek / first(Liczba_bialek), 2),
        Procent_peptydow_wzgledem_startu = round(100 * Liczba_peptydow / first(Liczba_peptydow), 2)
      ) %>%
      select(Wariant, Wariant_ID, Gradient, everything())
  }

  bind_rows(
    zrob("Cała FASTA", z$start_all, z$measured_all, z$rt_all, z$pep_all, z$frag_all),
    zrob("Panel 741 białek", z$start_panel, z$measured_panel, z$rt_panel, z$pep_panel, z$frag_panel)
  )
}

z_0_30 <- zbuduj_zbiory_gradient(E1_trypsyna_0DPMC, E2_0_30, panel_ids)
z_0_60 <- zbuduj_zbiory_gradient(E1_trypsyna_0DPMC, E2_0_60, panel_ids)
z_1_30 <- zbuduj_zbiory_gradient(E1_trypsyna_1DPMC, E2_30, panel_ids)
z_1_60 <- zbuduj_zbiory_gradient(E1_trypsyna_1DPMC, E2_60, panel_ids)

porownanie_gradienty_szczegolowe <- bind_rows(
  podsumuj_wariant_gradient(z_0_30, "Trypsyna 0 DPMC", "0DPMC", "30 min"),
  podsumuj_wariant_gradient(z_0_60, "Trypsyna 0 DPMC", "0DPMC", "60 min"),
  podsumuj_wariant_gradient(z_1_30, "Trypsyna 1 DPMC po ragged", "1DPMC", "30 min"),
  podsumuj_wariant_gradient(z_1_60, "Trypsyna 1 DPMC po ragged", "1DPMC", "60 min")
)

porownanie_gradienty_szerokie <- porownanie_gradienty_szczegolowe %>%
  select(Wariant, Wariant_ID, Gradient, Zbior, Etap_ID, Liczba_bialek, Liczba_peptydow) %>%
  pivot_wider(
    names_from = Etap_ID,
    values_from = c(Liczba_bialek, Liczba_peptydow),
    names_glue = "{.value}_etap_{Etap_ID}"
  ) %>%
  transmute(
    Wariant, Wariant_ID, Gradient, Zbior,
    Bialka_teoretyczne = Liczba_bialek_etap_1,
    Peptydy_teoretyczne = Liczba_peptydow_etap_1,
    Bialka_zmierzone = Liczba_bialek_etap_2,
    Peptydy_zmierzone = Liczba_peptydow_etap_2,
    Bialka_z_RT = Liczba_bialek_etap_3,
    Peptydy_z_RT = Liczba_peptydow_etap_3,
    Bialka_po_PEP = Liczba_bialek_etap_4,
    Peptydy_po_PEP = Liczba_peptydow_etap_4,
    Bialka_po_fragmentach = Liczba_bialek_etap_5,
    Peptydy_po_fragmentach = Liczba_peptydow_etap_5,
    Pokrycie_bialek_zmierzonych_proc = round(100 * Bialka_zmierzone / Bialka_teoretyczne, 2),
    Pokrycie_peptydow_zmierzonych_proc = round(100 * Peptydy_zmierzone / Peptydy_teoretyczne, 2),
    Pokrycie_bialek_po_fragmentach_proc = round(100 * Bialka_po_fragmentach / Bialka_teoretyczne, 2),
    Pokrycie_peptydow_po_fragmentach_proc = round(100 * Peptydy_po_fragmentach / Peptydy_teoretyczne, 2)
  )

porownanie_0_vs_1_roznice <- porownanie_gradienty_szerokie %>%
  select(
    Wariant_ID, Gradient, Zbior,
    Bialka_zmierzone, Peptydy_zmierzone,
    Bialka_po_PEP, Peptydy_po_PEP,
    Bialka_po_fragmentach, Peptydy_po_fragmentach
  ) %>%
  pivot_wider(
    names_from = Wariant_ID,
    values_from = c(
      Bialka_zmierzone, Peptydy_zmierzone,
      Bialka_po_PEP, Peptydy_po_PEP,
      Bialka_po_fragmentach, Peptydy_po_fragmentach
    ),
    names_sep = "__"
  ) %>%
  mutate(
    Roznica_bialek_zmierzonych_1_minus_0 = Bialka_zmierzone__1DPMC - Bialka_zmierzone__0DPMC,
    Roznica_peptydow_zmierzonych_1_minus_0 = Peptydy_zmierzone__1DPMC - Peptydy_zmierzone__0DPMC,
    Roznica_bialek_po_PEP_1_minus_0 = Bialka_po_PEP__1DPMC - Bialka_po_PEP__0DPMC,
    Roznica_peptydow_po_PEP_1_minus_0 = Peptydy_po_PEP__1DPMC - Peptydy_po_PEP__0DPMC,
    Roznica_bialek_po_fragmentach_1_minus_0 = Bialka_po_fragmentach__1DPMC - Bialka_po_fragmentach__0DPMC,
    Roznica_peptydow_po_fragmentach_1_minus_0 = Peptydy_po_fragmentach__1DPMC - Peptydy_po_fragmentach__0DPMC
  )

# Zbiory porównawcze nie są potrzebne w dalszej głównej selekcji.
rm(E2_0_30, E2_0_60, z_0_30, z_0_60, z_1_30, z_1_60); gc()

# ------------------------------------------------------------------------------
# 7. GŁÓWNA SELEKCJA — TRYPSYNA 1 DPMC PO RAGGED, DOCELOWY RT30
# ------------------------------------------------------------------------------

# Kolejność filtrów jest zachowana również w tabelach redukcji, dzięki czemu można śledzić wpływ każdego kryterium.

zbuduj_zbiory_glowne <- function(start_df, E2_30_df, E2_60_df, panel_ids) {
  start_all <- start_df %>% distinct(Bialko, Peptyd, .keep_all = TRUE)
  start_panel <- start_all %>% filter(Bialko %in% panel_ids)

  zmierzone_30_60_ids <- bind_rows(
    E2_30_df %>% filter(Zmierzony) %>% select(Bialko, Peptyd),
    E2_60_df %>% filter(Zmierzony) %>% select(Bialko, Peptyd)
  ) %>%
    distinct(Bialko, Peptyd)

  measured_any_all <- start_all %>%
    inner_join(zmierzone_30_60_ids, by = c("Bialko", "Peptyd"))
  measured_any_panel <- measured_any_all %>% filter(Bialko %in% panel_ids)

  rt30_all <- E2_30_df %>%
    filter(
      Zmierzony,
      !is.na(RT_sredni),
      Wybrany_ladunek %in% dozwolone_ladunki_prekursora
    ) %>%
    distinct(Bialko, Peptyd, .keep_all = TRUE)
  rt30_panel <- rt30_all %>% filter(Bialko %in% panel_ids)

  pep_all <- rt30_all %>%
    filter(!is.na(Najlepszy_PEP), Najlepszy_PEP <= prog_PEP)
  pep_panel <- pep_all %>% filter(Bialko %in% panel_ids)

  frag_all <- pep_all %>%
    filter(!is.na(Max_dobrych_fragmentow), Max_dobrych_fragmentow >= min_dobrych_fragmentow)
  frag_panel <- frag_all %>% filter(Bialko %in% panel_ids)

  top2_all <- ranking_top_n(frag_all, n = 2L)
  top2_panel <- ranking_top_n(frag_panel, n = 2L)
  top3_all <- ranking_top_n(frag_all, n = 3L)
  top3_panel <- ranking_top_n(frag_panel, n = 3L)

  list(
    start_all = start_all, start_panel = start_panel,
    measured_any_all = measured_any_all, measured_any_panel = measured_any_panel,
    rt30_all = rt30_all, rt30_panel = rt30_panel,
    pep_all = pep_all, pep_panel = pep_panel,
    frag_all = frag_all, frag_panel = frag_panel,
    top2_all = top2_all, top2_panel = top2_panel,
    top3_all = top3_all, top3_panel = top3_panel
  )
}

zbiory <- zbuduj_zbiory_glowne(E1_trypsyna_1DPMC, E2_30, E2_60, panel_ids)

E2_30_przed <- bind_rows(E2_30, E2_30_ragged_usuniete) %>%
  semi_join(klucze_E1_1_przed, by = c("Bialko", "Peptyd"))
E2_60_przed <- bind_rows(E2_60, E2_60_ragged_usuniete) %>%
  semi_join(klucze_E1_1_przed, by = c("Bialko", "Peptyd"))
zbiory_przed <- zbuduj_zbiory_glowne(
  E1_trypsyna_1DPMC_przed,
  E2_30_przed,
  E2_60_przed,
  panel_ids
)
rm(E2_30_przed, E2_60_przed, E2_30_ragged_usuniete, E2_60_ragged_usuniete); gc()

list2env(zbiory, envir = .GlobalEnv)

# Zestawia liczebności peptydów i białek po kolejnych kryteriach głównej selekcji.
zrob_redukcje_glowne <- function(z, etykieta_startu) {
  zrob <- function(zbior, start, measured, rt, pep, frag, top3) {
    bind_rows(
      licz_podsumowanie(start, zbior, 1L, etykieta_startu, "Teoretyczne peptydy unikatowe z E1 HUMAN-ONLY."),
      licz_podsumowanie(measured, zbior, 2L, "Zmierzone w DIA-NN 30/60 min", "Peptydy zmierzone w co najmniej jednym gradiencie."),
      licz_podsumowanie(rt, zbior, 3L, "Obserwowany RT30 i ładunek 2+/3+", "Peptydy z obserwowanym RT w docelowym gradiencie 30 min."),
      licz_podsumowanie(pep, zbior, 4L, "PEP <= 0,01", "PEP jest jedynym kryterium jakości identyfikacji."),
      licz_podsumowanie(frag, zbior, 5L, "PEP <= 0,01 i >=3 kwalifikujące się fragmenty", "Minimum trzy kwalifikujące się fragmenty na podstawie danych EXPORT-QUANT."),
      licz_podsumowanie(top3, zbior, 6L, "Top 3 peptydy na białko", "Maksymalnie trzy najlepsze peptydy na reprezentowane białko.")
    ) %>%
      mutate(
        Procent_bialek_wzgledem_startu = round(100 * Liczba_bialek / first(Liczba_bialek), 2),
        Procent_peptydow_wzgledem_startu = round(100 * Liczba_peptydow / first(Liczba_peptydow), 2)
      )
  }

  list(
    all = zrob("Cała FASTA", z$start_all, z$measured_any_all, z$rt30_all,
               z$pep_all, z$frag_all, z$top3_all),
    panel = zrob("Panel 741 białek", z$start_panel, z$measured_any_panel, z$rt30_panel,
                 z$pep_panel, z$frag_panel, z$top3_panel)
  )
}

redukcje <- zrob_redukcje_glowne(zbiory, "Start: Trypsyna 1 DPMC po ragged")
redukcja_all <- redukcje$all
redukcja_panel <- redukcje$panel
redukcja_laczna <- bind_rows(redukcja_all, redukcja_panel)

redukcje_przed <- zrob_redukcje_glowne(zbiory_przed, "Start: Trypsyna 1 DPMC przed ragged")
porownanie_przed_po_ragged <- bind_rows(
  bind_rows(redukcje_przed$all, redukcje_przed$panel) %>% mutate(Wersja = "Przed filtrem ragged"),
  redukcja_laczna %>% mutate(Wersja = "Po filtrze ragged")
) %>%
  select(Wersja, everything()) %>%
  group_by(Zbior, Etap_ID) %>%
  mutate(
    Zmiana_bialek_wobec_przed = Liczba_bialek - Liczba_bialek[Wersja == "Przed filtrem ragged"][1],
    Zmiana_peptydow_wobec_przed = Liczba_peptydow - Liczba_peptydow[Wersja == "Przed filtrem ragged"][1]
  ) %>%
  ungroup()

rm(zbiory_przed, E1_trypsyna_1DPMC_przed); gc()

porownanie_top2_top3 <- bind_rows(
  licz_podsumowanie(top2_all, "Cała FASTA", 2L, "Top 2", "Alternatywny limit maksymalnie dwóch peptydów na białko."),
  licz_podsumowanie(top3_all, "Cała FASTA", 3L, "Top 3", "Alternatywny limit maksymalnie trzech peptydów na białko."),
  licz_podsumowanie(top2_panel, "Panel 741 białek", 2L, "Top 2", "Alternatywny limit maksymalnie dwóch peptydów na białko."),
  licz_podsumowanie(top3_panel, "Panel 741 białek", 3L, "Top 3", "Alternatywny limit maksymalnie trzech peptydów na białko.")
) %>%
  rename(Limit_ID = Etap_ID, Limit = Etap)

# ------------------------------------------------------------------------------
# 8. DYSTRYBUCJE I STATUSY PANELU
# ------------------------------------------------------------------------------

dystr_all <- bind_rows(
  zrob_dystrybucje(start_all, wszystkie_bialka_FASTA_do_dystrybucji, "Cała FASTA", 1L, "Start po filtrze ragged"),
  zrob_dystrybucje(measured_any_all, wszystkie_bialka_FASTA_do_dystrybucji, "Cała FASTA", 2L, "Zmierzone w DIA-NN 30/60 min"),
  zrob_dystrybucje(rt30_all, wszystkie_bialka_FASTA_do_dystrybucji, "Cała FASTA", 3L, "Obserwowany RT30 i ładunek 2+/3+"),
  zrob_dystrybucje(pep_all, wszystkie_bialka_FASTA_do_dystrybucji, "Cała FASTA", 4L, "PEP <= 0,01"),
  zrob_dystrybucje(frag_all, wszystkie_bialka_FASTA_do_dystrybucji, "Cała FASTA", 5L, "PEP <= 0,01 i >=3 fragmenty"),
  zrob_dystrybucje(top3_all, wszystkie_bialka_FASTA_do_dystrybucji, "Cała FASTA", 6L, "Top 3 peptydy na białko")
)

dystr_panel <- bind_rows(
  zrob_dystrybucje(start_panel, panel_ids, "Panel 741 białek", 1L, "Start po filtrze ragged"),
  zrob_dystrybucje(measured_any_panel, panel_ids, "Panel 741 białek", 2L, "Zmierzone w DIA-NN 30/60 min"),
  zrob_dystrybucje(rt30_panel, panel_ids, "Panel 741 białek", 3L, "Obserwowany RT30 i ładunek 2+/3+"),
  zrob_dystrybucje(pep_panel, panel_ids, "Panel 741 białek", 4L, "PEP <= 0,01"),
  zrob_dystrybucje(frag_panel, panel_ids, "Panel 741 białek", 5L, "PEP <= 0,01 i >=3 fragmenty"),
  zrob_dystrybucje(top3_panel, panel_ids, "Panel 741 białek", 6L, "Top 3 peptydy na białko")
)
dystr_laczna <- bind_rows(dystr_all, dystr_panel)

liczby_na_bialko <- function(df, nazwa) {
  df %>%
    distinct(Bialko, Peptyd) %>%
    count(Bialko, name = nazwa)
}

status_panel <- tibble(Bialko = panel_ids) %>%
  left_join(liczby_na_bialko(start_panel, "N_start"), by = "Bialko") %>%
  left_join(liczby_na_bialko(measured_any_panel, "N_zmierzone_30_60"), by = "Bialko") %>%
  left_join(liczby_na_bialko(rt30_panel, "N_RT30"), by = "Bialko") %>%
  left_join(liczby_na_bialko(pep_panel, "N_PEP"), by = "Bialko") %>%
  left_join(liczby_na_bialko(frag_panel, "N_fragmenty"), by = "Bialko") %>%
  left_join(liczby_na_bialko(top2_panel, "N_top2"), by = "Bialko") %>%
  left_join(liczby_na_bialko(top3_panel, "N_top3"), by = "Bialko") %>%
  mutate(
    across(starts_with("N_"), ~ replace_na(.x, 0L)),
    Status_koncowy = case_when(
      N_top3 > 0 ~ "Reprezentowane na liście Top 3",
      N_fragmenty > 0 ~ "Spełnia kryteria, lecz błąd rankingu Top 3",
      N_PEP > 0 ~ "PEP <= 0,01, ale mniej niż 3 kwalifikujące się fragmenty EXPORT-QUANT",
      N_RT30 > 0 ~ "Obserwowany RT30, ale PEP > 0,01 lub brak PEP",
      N_zmierzone_30_60 > 0 ~ "Zmierzone w 30/60 min, ale brak prawidłowego RT30 lub ładunku 2+/3+",
      N_start > 0 ~ "Peptyd teoretyczny istnieje, ale brak detekcji DIA-NN 30/60",
      TRUE ~ "Brak unikatowego peptydu Trypsyna 1 DPMC po ragged"
    )
  )

status_panel_podsumowanie <- status_panel %>%
  count(Status_koncowy, name = "Liczba_bialek") %>%
  mutate(Procent_panelu = round(100 * Liczba_bialek / length(panel_ids), 2))

# ------------------------------------------------------------------------------
# 9. KONTROLE STRUKTURALNE — BEZ STARYCH WARTOŚCI WPISANYCH NA SZTYWNO
# ------------------------------------------------------------------------------

# Kontrole porównują strukturę i liczebności wygenerowanych danych z wartościami odczytanymi z bieżącego uruchomienia.

wiersz_E1_0 <- podsumowanie_E1 %>% filter(str_detect(Wariant, "Trypsyna, 0 DPMC")) %>% slice(1)
wiersz_E1_1 <- podsumowanie_E1 %>% filter(str_detect(Wariant, "Trypsyna, 1 DPMC po filtrze ragged")) %>% slice(1)

# Kontrola specyficzna dla v3: jeżeli istnieje co najmniej jeden prekursor
# spełniający PEP <= 0,01, wybrany prekursor również musi spełniać PEP.
kontrola_pep_first <- bind_rows(
  E2_30 %>% mutate(Gradient_kontrola = "30min"),
  E2_60 %>% mutate(Gradient_kontrola = "60min")
) %>%
  filter(!is.na(N_prekursorow_PEP_OK), N_prekursorow_PEP_OK > 0L) %>%
  mutate(
    Naruszenie = is.na(Najlepszy_PEP) | Najlepszy_PEP > prog_PEP
  )

# Dodaje pojedynczy wynik kontroli do tabeli testów poprawności.
kontrola <- function(nazwa, wynik, obserwowana, oczekiwana) {
  tibble(
    Kontrola = nazwa,
    Wynik = isTRUE(wynik),
    Wartosc_obserwowana = as.character(obserwowana),
    Wartosc_oczekiwana = as.character(oczekiwana)
  )
}

kontrole <- bind_rows(
  kontrola("Mianownik FASTA = 20 663 Homo sapiens",
           liczba_bialek_FASTA == oczekiwana_liczba_bialek_FASTA,
           liczba_bialek_FASTA, oczekiwana_liczba_bialek_FASTA),
  kontrola("Panel zawiera 741 unikatowych identyfikatorów",
           length(panel_ids) == oczekiwana_liczba_bialek_panelu,
           length(panel_ids), oczekiwana_liczba_bialek_panelu),
  kontrola("E1 0 DPMC zgadza się z podsumowaniem ETAPU 00/01 — peptydy",
           nrow(E1_trypsyna_0DPMC) == as.integer(na_liczbe(wiersz_E1_0$Peptydy_unikatowe)),
           nrow(E1_trypsyna_0DPMC), as.integer(na_liczbe(wiersz_E1_0$Peptydy_unikatowe))),
  kontrola("E1 1 DPMC zgadza się z podsumowaniem ETAPU 00/01 — peptydy",
           nrow(E1_trypsyna_1DPMC) == as.integer(na_liczbe(wiersz_E1_1$Peptydy_unikatowe)),
           nrow(E1_trypsyna_1DPMC), as.integer(na_liczbe(wiersz_E1_1$Peptydy_unikatowe))),
  kontrola("E1 0 DPMC zgadza się z podsumowaniem ETAPU 00/01 — białka",
           n_distinct(E1_trypsyna_0DPMC$Bialko) == as.integer(na_liczbe(wiersz_E1_0$Bialka_z_min_1_peptydem_unikatowym)),
           n_distinct(E1_trypsyna_0DPMC$Bialko), as.integer(na_liczbe(wiersz_E1_0$Bialka_z_min_1_peptydem_unikatowym))),
  kontrola("E1 1 DPMC zgadza się z podsumowaniem ETAPU 00/01 — białka",
           n_distinct(E1_trypsyna_1DPMC$Bialko) == as.integer(na_liczbe(wiersz_E1_1$Bialka_z_min_1_peptydem_unikatowym)),
           n_distinct(E1_trypsyna_1DPMC$Bialko), as.integer(na_liczbe(wiersz_E1_1$Bialka_z_min_1_peptydem_unikatowym))),
  kontrola("Panel E1 0 DPMC zgadza się z podsumowaniem — peptydy",
           nrow(filter(E1_trypsyna_0DPMC, Bialko %in% panel_ids)) == as.integer(na_liczbe(wiersz_E1_0$Panel_peptydy_unikatowe)),
           nrow(filter(E1_trypsyna_0DPMC, Bialko %in% panel_ids)), as.integer(na_liczbe(wiersz_E1_0$Panel_peptydy_unikatowe))),
  kontrola("Panel E1 1 DPMC zgadza się z podsumowaniem — peptydy",
           nrow(filter(E1_trypsyna_1DPMC, Bialko %in% panel_ids)) == as.integer(na_liczbe(wiersz_E1_1$Panel_peptydy_unikatowe)),
           nrow(filter(E1_trypsyna_1DPMC, Bialko %in% panel_ids)), as.integer(na_liczbe(wiersz_E1_1$Panel_peptydy_unikatowe))),
  kontrola("Brak technicznego rekordu Biognosys iRT w E1",
           !any(str_detect(E1_trypsyna_0DPMC$Bialko, wzorzec_standardu_iRT)) &&
             !any(str_detect(E1_trypsyna_1DPMC$Bialko, wzorzec_standardu_iRT)),
           sum(str_detect(c(E1_trypsyna_0DPMC$Bialko, E1_trypsyna_1DPMC$Bialko), wzorzec_standardu_iRT)), 0),
  kontrola("Brak ragged w finalnym E1 Trypsyna 1 DPMC",
           !any(czy_ragged(E1_trypsyna_1DPMC$Peptyd)),
           sum(czy_ragged(E1_trypsyna_1DPMC$Peptyd)), 0),
  kontrola("Brak metioniny w finalnym E1 Trypsyna 1 DPMC",
           !any(str_detect(E1_trypsyna_1DPMC$Peptyd, "M")),
           sum(str_detect(E1_trypsyna_1DPMC$Peptyd, "M")), 0),
  kontrola("Długości E1 Trypsyna 1 DPMC mieszczą się w 7-25 aa",
           all(nchar(E1_trypsyna_1DPMC$Peptyd) >= 7 & nchar(E1_trypsyna_1DPMC$Peptyd) <= 25),
           paste0(min(nchar(E1_trypsyna_1DPMC$Peptyd)), "-", max(nchar(E1_trypsyna_1DPMC$Peptyd))), "7-25"),
  kontrola("Globalna unikatowość E1 0 DPMC",
           kontrola_globalnej_unikatowosci(E1_trypsyna_0DPMC), "sprawdzono", "TRUE"),
  kontrola("Globalna unikatowość E1 1 DPMC",
           kontrola_globalnej_unikatowosci(E1_trypsyna_1DPMC), "sprawdzono", "TRUE"),
  kontrola("Każda tabela E2 pokrywa 100% aktualnego E1",
           all(kontrola_pokrycia_E2$Brakujace_pary_E1_w_E2 == 0L),
           sum(kontrola_pokrycia_E2$Brakujace_pary_E1_w_E2), 0),
  kontrola("Każdy raport EXPORT-QUANT ma sloty Fr.0-Fr.11",
           all(kontrola_export_quant_run$Min_slot == 0L & kontrola_export_quant_run$Max_slot == 11L),
           paste(kontrola_export_quant_run$Min_slot, kontrola_export_quant_run$Max_slot, sep = "-", collapse = "; "), "0-11 dla każdego runu"),
  kontrola("Parser Fr.Id działa dla wszystkich niepustych identyfikatorów fragmentów",
           all(kontrola_export_quant_run$N_parser_nieudany == 0L),
           sum(kontrola_export_quant_run$N_parser_nieudany), 0),
  kontrola("Brak NA w Fr.Quantity i Fr.Score",
           all(kontrola_export_quant_run$N_quantity_NA == 0L & kontrola_export_quant_run$N_score_NA == 0L),
           paste0("Quantity NA=", sum(kontrola_export_quant_run$N_quantity_NA), "; Score NA=", sum(kontrola_export_quant_run$N_score_NA)), "0 / 0"),
  kontrola("Top fragmenty są jonami y 1+/2+ bez adnotacji utraty i z m/z >300",
           all(tolower(mapa_top_fragmentow$Fragment.Type) %in% dozwolone_typy_fragmentow &
               mapa_top_fragmentow$Fragment.Charge %in% dozwolone_ladunki_fragmentow &
               mapa_top_fragmentow$Fragment.Annotation == "" &
               mapa_top_fragmentow$Product.Mz > minimalne_product_mz),
           sum(!(tolower(mapa_top_fragmentow$Fragment.Type) %in% dozwolone_typy_fragmentow &
                 mapa_top_fragmentow$Fragment.Charge %in% dozwolone_ladunki_fragmentow &
                 mapa_top_fragmentow$Fragment.Annotation == "" &
                 mapa_top_fragmentow$Product.Mz > minimalne_product_mz)), 0),
  kontrola("Wszystkie rekordy RT30 mają obserwowany RT",
           all(!is.na(rt30_all$RT_sredni)), sum(is.na(rt30_all$RT_sredni)), 0),
  kontrola("Wszystkie rekordy RT30 mają ładunek 2+ lub 3+",
           all(rt30_all$Wybrany_ladunek %in% dozwolone_ladunki_prekursora),
           sum(!rt30_all$Wybrany_ladunek %in% dozwolone_ladunki_prekursora), 0),
  kontrola("Wszystkie rekordy RT30 mają teoretyczne m/z prekursora do PRM",
           all(is.finite(rt30_all$Wybrane_mz)),
           sum(!is.finite(rt30_all$Wybrane_mz)), 0),
  kontrola("Wszystkie rekordy po PEP spełniają PEP <= 0,01",
           all(!is.na(pep_all$Najlepszy_PEP) & pep_all$Najlepszy_PEP <= prog_PEP),
           sum(is.na(pep_all$Najlepszy_PEP) | pep_all$Najlepszy_PEP > prog_PEP), 0),
  kontrola("PEP-first: fragmenty nie wpływają na przejście etapu PEP",
           !any(kontrola_pep_first$Naruszenie),
           sum(kontrola_pep_first$Naruszenie), 0),
  kontrola("Wszystkie rekordy po fragmentach mają >=3 kwalifikujące się fragmenty",
           all(!is.na(frag_all$Max_dobrych_fragmentow) & frag_all$Max_dobrych_fragmentow >= min_dobrych_fragmentow),
           sum(is.na(frag_all$Max_dobrych_fragmentow) | frag_all$Max_dobrych_fragmentow < min_dobrych_fragmentow), 0),
  kontrola("Top 2: maksymalnie 2 peptydy na białko",
           max(count(top2_all, Bialko)$n) <= 2L, max(count(top2_all, Bialko)$n), 2),
  kontrola("Top 3: maksymalnie 3 peptydy na białko",
           max(count(top3_all, Bialko)$n) <= 3L, max(count(top3_all, Bialko)$n), 3),
  kontrola("Top 2 jest podzbiorem Top 3",
           nrow(anti_join(top2_all, top3_all, by = c("Bialko", "Peptyd"))) == 0L,
           nrow(anti_join(top2_all, top3_all, by = c("Bialko", "Peptyd"))), 0),
  kontrola("Top 3 jest podzbiorem kandydatów po fragmentach",
           nrow(anti_join(top3_all, frag_all, by = c("Bialko", "Peptyd"))) == 0L,
           nrow(anti_join(top3_all, frag_all, by = c("Bialko", "Peptyd"))), 0),
  kontrola("Top 2 zachowuje wszystkie białka po kryterium fragmentów",
           n_distinct(top2_all$Bialko) == n_distinct(frag_all$Bialko),
           n_distinct(top2_all$Bialko), n_distinct(frag_all$Bialko)),
  kontrola("Top 3 zachowuje wszystkie białka po kryterium fragmentów",
           n_distinct(top3_all$Bialko) == n_distinct(frag_all$Bialko),
           n_distinct(top3_all$Bialko), n_distinct(frag_all$Bialko)),
  kontrola("Eksportowane tabele jakościowe nie zawierają Q.Value",
           !any(str_detect(colnames(usun_QValue(frag_all)), regex("QValue|Q[.]Value|Q_Value", ignore_case = TRUE))),
           sum(str_detect(colnames(usun_QValue(frag_all)), regex("QValue|Q[.]Value|Q_Value", ignore_case = TRUE))), 0),
  kontrola("Redukcja peptydów jest monotoniczna — FASTA",
           all(diff(redukcja_all$Liczba_peptydow) <= 0),
           all(diff(redukcja_all$Liczba_peptydow) <= 0), TRUE),
  kontrola("Redukcja białek jest monotoniczna — FASTA",
           all(diff(redukcja_all$Liczba_bialek) <= 0),
           all(diff(redukcja_all$Liczba_bialek) <= 0), TRUE),
  kontrola("Redukcja peptydów jest monotoniczna — panel",
           all(diff(redukcja_panel$Liczba_peptydow) <= 0),
           all(diff(redukcja_panel$Liczba_peptydow) <= 0), TRUE),
  kontrola("Każde z 741 białek ma dokładnie jeden status końcowy",
           nrow(status_panel) == 741L && n_distinct(status_panel$Bialko) == 741L && !any(is.na(status_panel$Status_koncowy)),
           paste(nrow(status_panel), n_distinct(status_panel$Bialko), sep = "/"), "741/741"),
  kontrola("Nie występuje niemożliwy status błędu rankingu Top 3",
           !any(status_panel$Status_koncowy == "Spełnia kryteria, lecz błąd rankingu Top 3"),
           sum(status_panel$Status_koncowy == "Spełnia kryteria, lecz błąd rankingu Top 3"), 0)
)

# Snapshot liczb jest zapisywany dynamicznie. Nie jest używany jako twardo wpisany
# warunek zatrzymujący skrypt.
snapshot_wynikow <- bind_rows(
  redukcja_laczna %>% mutate(Analiza = "Główna selekcja 1 DPMC"),
  porownanie_gradienty_szczegolowe %>%
    transmute(
      Zbior, Etap_ID, Etap, Opis,
      Liczba_bialek, Liczba_peptydow,
      Procent_bialek_wzgledem_startu,
      Procent_peptydow_wzgledem_startu,
      Analiza = paste(Wariant, Gradient)
    )
)

openxlsx::write.xlsx(
  list(
    Kontrole_strukturalne = kontrole,
    Export_quant_runy = kontrola_export_quant_run,
    Metoda_fragmentow = kontrola_fragmentow_metoda,
    Pokrycie_E2 = kontrola_pokrycia_E2,
    Statusy_panelu = status_panel_podsumowanie,
    Snapshot_wynikow = snapshot_wynikow
  ),
  file.path(folder_kontrola, "KONTROLA_ETAP_02_FINAL.xlsx"),
  overwrite = TRUE
)
zapisz_csv2(kontrole, "KONTROLA_ETAP_02_FINAL.csv")
zapisz_csv2(kontrola_pokrycia_E2, "KONTROLA_pokrycia_E2_aktualnym_E1.csv")
readr::write_excel_csv2(kontrola_export_quant_run, file.path(folder_kontrola, "KONTROLA_export_quant_runy.csv"))
readr::write_excel_csv2(kontrola_fragmentow_metoda, file.path(folder_kontrola, "KONTROLA_metoda_fragmentow_export_quant.csv"))

if (!all(kontrole$Wynik)) {
  stop(
    "Co najmniej jedna kontrola ETAPU 02 ma wynik FALSE. ",
    "Sprawdź 03_Kontrola/KONTROLA_ETAP_02_FINAL.xlsx."
  )
}

# ------------------------------------------------------------------------------
# 10. EKSPORT TABEL
# ------------------------------------------------------------------------------

metodologia_porownania <- tibble::tribble(
  ~Element, ~Opis,
  "Warianty", "Trypsyna 0 DPMC oraz Trypsyna 1 DPMC po filtrze ragged.",
  "Gradienty", "Dane E2 z gradientów 30 i 60 min analizowane oddzielnie.",
  "Proteom", "Wyłącznie 20 663 rekordy Homo sapiens; standard Biognosys iRT wykluczony.",
  "Jakość identyfikacji", "Wyłącznie PEP <= 0,01; Q.Value nie jest kryterium selekcji. PEP jest stosowane przed rankingiem fragmentacyjnym stanu naładowania.",
  "Fragmenty", "Dopiero wśród prekursorów spełniających PEP <= 0,01 wybierany jest najlepszy stan naładowania na podstawie EXPORT-QUANT; wymagane są minimum trzy kwalifikujące się fragmenty: y, 1+/2+, m/z >300, Quantity >0 w co najmniej jednym runie; Top 3 wybierane wg średniego Fr.Score bez arbitralnego progu Score.",
  "Dalsza lista PRM", "Główna lista podstawowa opiera się na obserwowanym RT30 dla Trypsyny 1 DPMC po ragged."
)

openxlsx::write.xlsx(
  list(
    Metodologia = metodologia_porownania,
    Podsumowanie_szerokie = porownanie_gradienty_szerokie,
    Etapy_szczegolowe = porownanie_gradienty_szczegolowe,
    Roznice_1_minus_0 = porownanie_0_vs_1_roznice
  ),
  file.path(folder_tabele, "Tabela_00_Porownanie_Trypsyna_0_vs_1_DPMC_30_60min.xlsx"),
  overwrite = TRUE
)

propozycja_kolejnosci <- tibble::tribble(
  ~Kolejnosc, ~Kryterium, ~Sposob_zastosowania, ~Uzasadnienie,
  1, "Trypsyna 1 DPMC po filtrze ragged", "E1 HUMAN-ONLY, długość 7-25 aa, bez metioniny, globalnie unikatowe sekwencje.", "Skorygowany teoretyczny zbiór startowy.",
  2, "Obecność w DIA-NN 30 lub 60 min", "Peptyd zmierzony w co najmniej jednym z gradientów.", "Potwierdzenie eksperymentalnej obserwowalności.",
  3, "Obserwowany RT30 i ładunek 2+/3+", "Do listy podstawowej przechodzą rekordy z obserwowanym RT w gradiencie docelowym.", "Nie zastępujemy obserwowanego RT30 wartością przewidywaną.",
  4, "PEP <= 0,01", "PEP jest jedynym kryterium jakości identyfikacji.", "Zgodność z uwagą Pani Profesor; Q.Value nie jest używane.",
  5, ">=3 kwalifikujące się fragmenty", "Po zastosowaniu PEP wybierany jest najlepszy stan naładowania; musi posiadać minimum trzy kwalifikujące się fragmenty EXPORT-QUANT.", "Fragmenty są oceniane na podstawie rzeczywistych Quantity i Fr.Score z głównego raportu DIA-NN.",
  6, "Ranking jakościowy", "Najpierw PEP <= 0,01; następnie wśród spełniających PEP prekursorów: liczba kwalifikujących się fragmentów -> średni Fr.Score Top 3 -> Quantity -> PEP jako rozstrzygnięcie dalszego remisu -> intensywność prekursora.", "Zapewnia kolejność PEP -> fragmenty i uniemożliwia metrykom fragmentacyjnym wpływanie na przejście etapu PEP.",
  7, "Alternatywny limit Top 2 lub Top 3", "Generowane są obie listy; Top 3 pozostaje wejściem do ETAPU 03.", "Pozwala odpowiedzieć na oba warianty rozważane przez Panią Profesor."
)

wb <- openxlsx::createWorkbook()
for (nazwa in c(
  "Propozycja_kolejnosci", "Redukcja_cala_FASTA", "Redukcja_panel",
  "Dystrybucja_cala_FASTA", "Dystrybucja_panel", "Przed_po_ragged",
  "Porownanie_Top2_Top3", "Statusy_panelu", "Kontrole"
)) openxlsx::addWorksheet(wb, nazwa)

openxlsx::writeData(wb, "Propozycja_kolejnosci", propozycja_kolejnosci)
openxlsx::writeData(wb, "Redukcja_cala_FASTA", redukcja_all)
openxlsx::writeData(wb, "Redukcja_panel", redukcja_panel)
openxlsx::writeData(wb, "Dystrybucja_cala_FASTA", dystr_all)
openxlsx::writeData(wb, "Dystrybucja_panel", dystr_panel)
openxlsx::writeData(wb, "Przed_po_ragged", porownanie_przed_po_ragged)
openxlsx::writeData(wb, "Porownanie_Top2_Top3", porownanie_top2_top3)
openxlsx::writeData(wb, "Statusy_panelu", status_panel_podsumowanie)
openxlsx::writeData(wb, "Kontrole", kontrole)
openxlsx::saveWorkbook(
  wb,
  file.path(folder_tabele, "Tabela_01_Redukcja_kryteriow_Trypsyna_1MC_PEP_podsumowanie.xlsx"),
  overwrite = TRUE
)

# Nazwy Tabel 06 i 07 są zachowane dla zgodności z ETAPEM 03.
zapisz_excel(usun_QValue(pep_all) %>% arrange(Bialko, Peptyd),
              file.path(folder_tabele, "Tabela_02_CALA_FASTA_kandydaci_po_PEP_0_01.xlsx"))
zapisz_excel(usun_QValue(pep_panel) %>% arrange(Bialko, Peptyd),
              file.path(folder_tabele, "Tabela_03_PANEL_kandydaci_po_PEP_0_01.xlsx"))
zapisz_excel(usun_QValue(frag_all) %>% arrange(Bialko, Peptyd),
              file.path(folder_tabele, "Tabela_04_CALA_FASTA_kandydaci_po_PEP_0_01_i_min_3_fragmentach.xlsx"))
zapisz_excel(usun_QValue(frag_panel) %>% arrange(Bialko, Peptyd),
              file.path(folder_tabele, "Tabela_05_PANEL_kandydaci_po_PEP_0_01_i_min_3_fragmentach.xlsx"))
zapisz_excel(usun_QValue(top3_all) %>% arrange(Bialko, Pozycja_w_bialku, Peptyd),
              file.path(folder_tabele, "Tabela_06_CALA_FASTA_top3_peptydy_na_bialko_PEP.xlsx"))
zapisz_excel(usun_QValue(top3_panel) %>% arrange(Bialko, Pozycja_w_bialku, Peptyd),
              file.path(folder_tabele, "Tabela_07_PANEL_top3_peptydy_na_bialko_PEP.xlsx"))
zapisz_excel(usun_QValue(top2_all) %>% arrange(Bialko, Pozycja_w_bialku, Peptyd),
              file.path(folder_tabele, "Tabela_08_CALA_FASTA_top2_peptydy_na_bialko_PEP.xlsx"))
zapisz_excel(usun_QValue(top2_panel) %>% arrange(Bialko, Pozycja_w_bialku, Peptyd),
              file.path(folder_tabele, "Tabela_09_PANEL_top2_peptydy_na_bialko_PEP.xlsx"))
zapisz_excel(status_panel,
              file.path(folder_kontrola, "Kontrola_PANEL_status_bialek_po_kryteriach_PEP.xlsx"))

zapisz_csv2(porownanie_gradienty_szerokie, "ETAP_02_porownanie_0_vs_1_DPMC_30_60_szerokie.csv")
zapisz_csv2(porownanie_gradienty_szczegolowe, "ETAP_02_porownanie_0_vs_1_DPMC_30_60_szczegolowe.csv")
zapisz_csv2(porownanie_0_vs_1_roznice, "ETAP_02_roznice_1DPMC_minus_0DPMC.csv")
zapisz_csv2(redukcja_all, "ETAP_02_redukcja_CALA_FASTA.csv")
zapisz_csv2(redukcja_panel, "ETAP_02_redukcja_PANEL_741.csv")
zapisz_csv2(dystr_all, "ETAP_02_dystrybucja_CALA_FASTA.csv")
zapisz_csv2(dystr_panel, "ETAP_02_dystrybucja_PANEL_741.csv")
zapisz_csv2(porownanie_przed_po_ragged, "ETAP_02_porownanie_przed_po_ragged.csv")
zapisz_csv2(porownanie_top2_top3, "ETAP_02_porownanie_Top2_Top3.csv")
zapisz_csv2(status_panel, "ETAP_02_status_741_bialek_panelu.csv")

# ------------------------------------------------------------------------------
# 11. WYKRESY
# ------------------------------------------------------------------------------

porownanie_plot <- porownanie_gradienty_szczegolowe %>%
  filter(Etap_ID %in% c(2L, 5L)) %>%
  mutate(
    Etap_skrocony = recode(
      Etap,
      "Zmierzone w danym gradiencie" = "Zmierzone",
      "PEP <= 0,01 i >=3 fragmenty" = "PEP i >=3 fragmenty"
    ),
    Grupa = paste(Gradient, Etap_skrocony, sep = " — ")
  )

p_por_pep <- ggplot(
  porownanie_plot,
  aes(x = Grupa, y = Liczba_peptydow, fill = Wariant)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(
    aes(label = scales::label_number(big.mark = " ")(Liczba_peptydow)),
    position = position_dodge(width = 0.8), vjust = -0.2, size = 3
  ) +
  facet_wrap(~Zbior, scales = "free_y") +
  scale_y_continuous(labels = scales::label_number(big.mark = " "),
                     expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "Porównanie liczby peptydów: Trypsyna 0 i 1 DPMC",
    subtitle = "Gradienty 30 i 60 min; 1 DPMC po filtrze ragged; proteom HUMAN-ONLY",
    x = "Gradient i etap",
    y = "Liczba peptydów",
    fill = "Wariant trawienia"
  ) +
  theme_prm(base_size = 11)
zapisz_wykres(p_por_pep, "Wykres_00_Porownanie_0_vs_1_DPMC_30_60_peptydy", 14, 8)

p_por_bialka <- ggplot(
  porownanie_plot,
  aes(x = Grupa, y = Liczba_bialek, fill = Wariant)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(
    aes(label = scales::label_number(big.mark = " ")(Liczba_bialek)),
    position = position_dodge(width = 0.8), vjust = -0.2, size = 3
  ) +
  facet_wrap(~Zbior, scales = "free_y") +
  scale_y_continuous(labels = scales::label_number(big.mark = " "),
                     expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "Porównanie liczby reprezentowanych białek: Trypsyna 0 i 1 DPMC",
    subtitle = "Gradienty 30 i 60 min; 1 DPMC po filtrze ragged; proteom HUMAN-ONLY",
    x = "Gradient i etap",
    y = "Liczba białek",
    fill = "Wariant trawienia"
  ) +
  theme_prm(base_size = 11)
zapisz_wykres(p_por_bialka, "Wykres_00B_Porownanie_0_vs_1_DPMC_30_60_bialka", 14, 8)

rysuj_redukcje <- function(df, zbior, metryka, nazwa_pliku) {
  d <- df %>%
    filter(Zbior == zbior) %>%
    mutate(
      Etap = forcats::fct_reorder(Etap, Etap_ID),
      Wartosc = .data[[metryka]],
      Etykieta = scales::label_number(big.mark = " ")(Wartosc)
    )

  p <- ggplot(d, aes(x = Etap, y = Wartosc)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = Etykieta), vjust = -0.25, size = 3.4) +
    scale_y_continuous(labels = scales::label_number(big.mark = " "),
                       expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = paste0("Zmiana ", ifelse(metryka == "Liczba_peptydow", "liczby peptydów", "liczby białek"), " po kryteriach"),
      subtitle = paste0(zbior, "; Trypsyna 1 DPMC po ragged; docelowy RT30; kryterium identyfikacji: PEP"),
      x = "Etap analizy",
      y = ifelse(metryka == "Liczba_peptydow", "Liczba peptydów", "Liczba białek")
    ) +
    theme_prm()
  zapisz_wykres(p, nazwa_pliku, 11, 6.2)
}

rysuj_redukcje(redukcja_laczna, "Cała FASTA", "Liczba_peptydow", "Wykres_01_CALA_FASTA_redukcja_liczby_peptydow_PEP")
rysuj_redukcje(redukcja_laczna, "Cała FASTA", "Liczba_bialek", "Wykres_02_CALA_FASTA_redukcja_liczby_bialek_PEP")
rysuj_redukcje(redukcja_laczna, "Panel 741 białek", "Liczba_peptydow", "Wykres_03_PANEL_redukcja_liczby_peptydow_PEP")
rysuj_redukcje(redukcja_laczna, "Panel 741 białek", "Liczba_bialek", "Wykres_04_PANEL_redukcja_liczby_bialek_PEP")

redukcja_long <- redukcja_laczna %>%
  select(Zbior, Etap_ID, Etap, Liczba_bialek, Liczba_peptydow) %>%
  pivot_longer(c(Liczba_bialek, Liczba_peptydow), names_to = "Metryka", values_to = "Wartosc") %>%
  mutate(
    Metryka = recode(Metryka, Liczba_bialek = "Białka", Liczba_peptydow = "Peptydy"),
    Etap = forcats::fct_reorder(Etap, Etap_ID),
    Etykieta = scales::label_number(big.mark = " ")(Wartosc)
  )

p_redukcja <- ggplot(redukcja_long, aes(x = Etap, y = Wartosc)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = Etykieta), vjust = -0.25, size = 3) +
  facet_grid(Metryka ~ Zbior, scales = "free_y") +
  scale_y_continuous(labels = scales::label_number(big.mark = " "),
                     expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "Redukcja liczby peptydów i białek po kryteriach eksperymentalnych",
    subtitle = "Trypsyna 1 DPMC po ragged; proteom HUMAN-ONLY; PEP <= 0,01",
    x = "Etap analizy", y = "Liczba obiektów"
  ) +
  theme_prm(base_size = 11)
zapisz_wykres(p_redukcja, "Wykres_05_Redukcja_laczna_FASTA_i_PANEL_PEP", 14, 9)

rysuj_dystrybucje <- function(df, zbior, nazwa_pliku) {
  d <- df %>%
    filter(Zbior == zbior) %>%
    mutate(
      Etap = forcats::fct_reorder(Etap, Etap_ID),
      Etykieta = paste0(scales::label_number(big.mark = " ")(Liczba_bialek), "\n", Procent, "%")
    )

  p <- ggplot(d, aes(x = Kategoria, y = Liczba_bialek)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = Etykieta), vjust = -0.25, size = 2.8) +
    facet_wrap(~Etap, scales = "free_y") +
    scale_y_continuous(labels = scales::label_number(big.mark = " "),
                       expand = expansion(mult = c(0, 0.22))) +
    labs(
      title = paste0("Dystrybucja liczby kandydackich peptydów na białko — ", zbior),
      subtitle = "Trypsyna 1 DPMC po ragged; kolejne kryteria eksperymentalne",
      x = "Liczba peptydów na białko", y = "Liczba białek"
    ) +
    theme_prm(base_size = 11)
  zapisz_wykres(p, nazwa_pliku, 13, 8)
}

rysuj_dystrybucje(dystr_laczna, "Cała FASTA", "Wykres_06_CALA_FASTA_dystrybucja_peptydow_po_kryteriach_PEP")
rysuj_dystrybucje(dystr_laczna, "Panel 741 białek", "Wykres_07_PANEL_dystrybucja_peptydow_po_kryteriach_PEP")

# ------------------------------------------------------------------------------
# 12. INDEKS I RAPORT
# ------------------------------------------------------------------------------

indeks_plikow <- tibble::tribble(
  ~Folder, ~Plik, ~Opis,
  "01_Tabele_Excel", "Tabela_00_Porownanie_Trypsyna_0_vs_1_DPMC_30_60min.xlsx", "Porównanie 0/1 DPMC dla 30/60 min — cała FASTA i panel.",
  "01_Tabele_Excel", "Tabela_01_Redukcja_kryteriow_Trypsyna_1MC_PEP_podsumowanie.xlsx", "Główne podsumowanie selekcji Trypsyna 1 DPMC po ragged.",
  "01_Tabele_Excel", "Tabela_02_CALA_FASTA_kandydaci_po_PEP_0_01.xlsx", "Kandydaci po PEP dla całej FASTA.",
  "01_Tabele_Excel", "Tabela_03_PANEL_kandydaci_po_PEP_0_01.xlsx", "Kandydaci po PEP dla panelu.",
  "01_Tabele_Excel", "Tabela_04_CALA_FASTA_kandydaci_po_PEP_0_01_i_min_3_fragmentach.xlsx", "Kandydaci po PEP i fragmentach — FASTA.",
  "01_Tabele_Excel", "Tabela_05_PANEL_kandydaci_po_PEP_0_01_i_min_3_fragmentach.xlsx", "Kandydaci po PEP i fragmentach — panel.",
  "01_Tabele_Excel", "Tabela_06_CALA_FASTA_top3_peptydy_na_bialko_PEP.xlsx", "Top 3 dla całej FASTA; wejście zgodne z ETAPEM 03.",
  "01_Tabele_Excel", "Tabela_07_PANEL_top3_peptydy_na_bialko_PEP.xlsx", "Top 3 dla panelu; wejście zgodne z ETAPEM 03.",
  "01_Tabele_Excel", "Tabela_08_CALA_FASTA_top2_peptydy_na_bialko_PEP.xlsx", "Alternatywna lista Top 2 dla całej FASTA.",
  "01_Tabele_Excel", "Tabela_09_PANEL_top2_peptydy_na_bialko_PEP.xlsx", "Alternatywna lista Top 2 dla panelu.",
  "02_Wykresy", "Wykres_00_Porownanie_0_vs_1_DPMC_30_60_peptydy.png/pdf", "Porównanie liczby peptydów 0/1 DPMC.",
  "02_Wykresy", "Wykres_00B_Porownanie_0_vs_1_DPMC_30_60_bialka.png/pdf", "Porównanie liczby białek 0/1 DPMC.",
  "03_Kontrola", "KONTROLA_ETAP_02_FINAL.xlsx", "Kontrole strukturalne, pokrycie E2, statusy panelu i dynamiczny snapshot liczb.",
  "03_Kontrola", "KONTROLA_export_quant_runy.csv", "Kontrola obecności i poprawności Fr.0-Fr.11 oraz parsera Fr.Id w czterech raportach.",
  "03_Kontrola", "KONTROLA_metoda_fragmentow_export_quant.csv", "Jawny zapis kryteriów kwalifikacji i rankingu fragmentów EXPORT-QUANT.",
  "03_Kontrola", "Kontrola_PANEL_status_bialek_po_kryteriach_PEP.xlsx", "Dokładny status wszystkich 741 białek panelowych.",
  "04_CSV", "ETAP_02_porownanie_0_vs_1_DPMC_30_60_szerokie.csv", "Zwięzłe porównanie wariantów i gradientów.",
  "04_CSV", "ETAP_02_porownanie_przed_po_ragged.csv", "Wpływ filtra ragged na główną selekcję."
)

zapisz_excel(
  indeks_plikow,
  file.path(folder_wyniki, "00_INDEKS_PLIKOW_CZESC_02_KRYTERIA_DOSWIADCZALNE_FINAL.xlsx")
)

raport <- capture.output({
  cat("ETAP 02 — KRYTERIA DOŚWIADCZALNE — EXPORT-QUANT — FINAL HUMAN-ONLY\n")
  cat("Data wykonania:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  cat("ETAP 00/01:", folder_etap01, "\n")
  cat("Źródło E2:", folder_dane_posrednie, "\n")
  cat("Źródło jakości fragmentów: główne raporty DIA-NN z --export-quant\n")
  cat("Panel:", plik_panel, "\n\n")
  cat("PORÓWNANIE 0/1 DPMC — PODSUMOWANIE\n")
  print(porownanie_gradienty_szerokie)
  cat("\nRÓŻNICE 1 DPMC MINUS 0 DPMC\n")
  print(porownanie_0_vs_1_roznice)
  cat("\nGŁÓWNA REDUKCJA — CAŁA FASTA\n")
  print(redukcja_all)
  cat("\nGŁÓWNA REDUKCJA — PANEL\n")
  print(redukcja_panel)
  cat("\nTOP 2 VS TOP 3\n")
  print(porownanie_top2_top3)
  cat("\nSTATUSY PANELU\n")
  print(status_panel_podsumowanie)
  cat("\nKONTROLE\n")
  print(kontrole)
})
writeLines(raport, file.path(folder_kontrola, "RAPORT_ETAP_02_FINAL.txt"), useBytes = TRUE)

# ------------------------------------------------------------------------------
# 13. PODSUMOWANIE W KONSOLI
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("GOTOWE — ETAP 02 EXPORT-QUANT FINAL HUMAN-ONLY\n")
cat("Folder wynikowy:\n", folder_wyniki, "\n")
cat("============================================================\n\n")
cat("Główna redukcja — cała FASTA:\n")
print(redukcja_all)
cat("\nGłówna redukcja — panel 741:\n")
print(redukcja_panel)
cat("\nPorównanie Top 2 i Top 3:\n")
print(porownanie_top2_top3)
cat("\nWszystkie kontrole: TRUE\n")
cat("Najważniejszy raport: 03_Kontrola/RAPORT_ETAP_02_FINAL.txt\n")

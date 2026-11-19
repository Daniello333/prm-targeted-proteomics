# ==============================================================================
# 03. PRZYGOTOWANIE LIST PEPTYDOWO-PREKURSOROWYCH
#
# Praca magisterska:
# „Planowanie celowanych eksperymentów proteomicznych z wykorzystaniem środowiska R”
# Autor: Daniel Drabik
#
# Skrypt porządkuje wyniki etapu 02 i przygotowuje listy peptydowo-prekursorowe
# Top 2 oraz Top 3 dla całej FASTA i panelu 741 białek. Nie wprowadza nowych
# progów jakościowych — zachowuje kryteria zastosowane wcześniej i sprawdza
# zgodność liczebności, strukturę tabel oraz brak peptydów ragged.
#
# Listy utworzone w tym etapie są wykorzystywane w analizie zatłoczenia czasu
# retencji oraz jako listy bazowe przy budowie wariantów rozszerzonych.
# ==============================================================================

pakiety <- c("tidyverse", "readxl", "openxlsx", "ggplot2", "scales")
for (p in pakiety) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, dependencies = TRUE)
  library(p, character.only = TRUE)
}
options(scipen = 999)

# ------------------------------------------------------------------------------
# 1. USTAWIENIA
# ------------------------------------------------------------------------------

folder_etap02 <- NULL

wzorzec_ragged_N <- "^(KK|KR|RK|RR)"
wzorzec_ragged_C <- "(KK|KR|RK|RR)$"

# Ustala katalog skryptu, w którym zostanie utworzony nowy folder wynikowy.
pobierz_katalog_skryptu <- function() {
  ramki <- sys.frames()
  for (i in rev(seq_along(ramki))) {
    ofile <- ramki[[i]]$ofile
    if (!is.null(ofile)) {
      return(dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE)))
    }
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

katalog_skryptu <- pobierz_katalog_skryptu()
znacznik <- format(Sys.time(), "%Y%m%d_%H%M%S")

folder_wyniki <- file.path(
  katalog_skryptu,
  paste0("WYNIKI_ETAP_03_LISTY_PEPTYDOWO_PREKURSOROWE_EXPORT_QUANT_HUMAN_ONLY_", znacznik)
)
folder_tabele <- file.path(folder_wyniki, "01_Tabele_Excel")
folder_csv <- file.path(folder_wyniki, "02_CSV_do_dalszego_przygotowania")
folder_wykresy <- file.path(folder_wyniki, "03_Wykresy")
folder_kontrola <- file.path(folder_wyniki, "04_Kontrola")
for (f in c(folder_wyniki, folder_tabele, folder_csv, folder_wykresy, folder_kontrola)) {
  dir.create(f, recursive = TRUE, showWarnings = FALSE)
}

korzenie_szukania <- unique(c(
  normalizePath(getwd(), winslash = "/", mustWork = FALSE),
  katalog_skryptu,
  normalizePath(path.expand("~/Downloads"), winslash = "/", mustWork = FALSE),
  normalizePath(path.expand("~/Desktop"), winslash = "/", mustWork = FALSE)
))
korzenie_szukania <- korzenie_szukania[dir.exists(korzenie_szukania)]

# ------------------------------------------------------------------------------
# 2. WYSZUKIWANIE ETAPU 02
# ------------------------------------------------------------------------------

wymagane_etap02 <- c(
  file.path("01_Tabele_Excel", "Tabela_06_CALA_FASTA_top3_peptydy_na_bialko_PEP.xlsx"),
  file.path("01_Tabele_Excel", "Tabela_07_PANEL_top3_peptydy_na_bialko_PEP.xlsx"),
  file.path("01_Tabele_Excel", "Tabela_08_CALA_FASTA_top2_peptydy_na_bialko_PEP.xlsx"),
  file.path("01_Tabele_Excel", "Tabela_09_PANEL_top2_peptydy_na_bialko_PEP.xlsx"),
  file.path("03_Kontrola", "Kontrola_PANEL_status_bialek_po_kryteriach_PEP.xlsx"),
  file.path("04_CSV", "ETAP_02_porownanie_Top2_Top3.csv")
)

# Sprawdza komplet plików wymaganych z etapu 02.
ma_pliki <- function(folder, wymagane) {
  dir.exists(folder) && all(file.exists(file.path(folder, wymagane)))
}

# Wybiera najnowszy spośród kilku zgodnych katalogów lub archiwów.
wybierz_najnowszy <- function(sciezki) {
  sciezki <- unique(sciezki[file.exists(sciezki) | dir.exists(sciezki)])
  if (length(sciezki) == 0) return(NULL)
  sciezki[which.max(file.info(sciezki)$mtime)]
}

lista_folderow <- function() {
  unique(unlist(lapply(korzenie_szukania, function(root) {
    c(root, list.dirs(root, recursive = TRUE, full.names = TRUE))
  }), use.names = FALSE))
}

lista_zipow <- function() {
  unique(unlist(lapply(korzenie_szukania, function(root) {
    list.files(root, pattern = "[.]zip$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  }), use.names = FALSE))
}

zip_zawiera <- function(zip_path, wymagane) {
  listing <- tryCatch(utils::unzip(zip_path, list = TRUE), error = function(e) NULL)
  if (is.null(listing)) return(FALSE)
  nazwy <- gsub("\\\\", "/", listing$Name)
  all(vapply(wymagane, function(rel) {
    any(endsWith(nazwy, gsub("\\\\", "/", rel)))
  }, logical(1)))
}

wypakuj_zip <- function(zip_path, wymagane) {
  tmp <- tempfile("ETAP03_ETAP02_")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  utils::unzip(zip_path, exdir = tmp)
  kandydaci <- c(tmp, list.dirs(tmp, recursive = TRUE, full.names = TRUE))
  kandydaci <- kandydaci[vapply(kandydaci, ma_pliki, logical(1), wymagane = wymagane)]
  wybierz_najnowszy(kandydaci)
}

# Odnajduje kompletny wynik etapu 02 albo wypakowuje go z archiwum ZIP.
znajdz_etap02 <- function() {
  foldery <- lista_folderow()
  kandydaci <- foldery[vapply(foldery, ma_pliki, logical(1), wymagane = wymagane_etap02)]

  preferowane <- kandydaci[stringr::str_detect(
    basename(kandydaci),
    stringr::regex(
      "^WYNIKI_ETAP_02_KRYTERIA_DOSWIADCZALNE_PO_RAGGED_EXPORT_QUANT_HUMAN_ONLY_",
      ignore_case = TRUE
    )
  )]
  if (length(preferowane) > 0) return(wybierz_najnowszy(preferowane))
  if (length(kandydaci) > 0) return(wybierz_najnowszy(kandydaci))

  zipy <- lista_zipow()
  zipy <- zipy[vapply(zipy, zip_zawiera, logical(1), wymagane = wymagane_etap02)]
  preferowane_zip <- zipy[stringr::str_detect(
    basename(zipy),
    stringr::regex("EXPORT_QUANT_HUMAN_ONLY", ignore_case = TRUE)
  )]
  if (length(preferowane_zip) > 0) zipy <- preferowane_zip
  if (length(zipy) == 0) return(NULL)

  wypakuj_zip(wybierz_najnowszy(zipy), wymagane_etap02)
}

if (is.null(folder_etap02)) folder_etap02 <- znajdz_etap02()
if (is.null(folder_etap02)) {
  stop("Nie znaleziono finalnego ETAPU 02 EXPORT-QUANT. Ustaw folder_etap02 ręcznie.")
}
folder_etap02 <- normalizePath(folder_etap02, winslash = "/", mustWork = TRUE)

cat("ETAP 02:\n", folder_etap02, "\n\n")
cat("Wyniki ETAPU 03:\n", folder_wyniki, "\n\n")

sciezki <- list(
  FASTA_Top2 = file.path(folder_etap02, "01_Tabele_Excel", "Tabela_08_CALA_FASTA_top2_peptydy_na_bialko_PEP.xlsx"),
  FASTA_Top3 = file.path(folder_etap02, "01_Tabele_Excel", "Tabela_06_CALA_FASTA_top3_peptydy_na_bialko_PEP.xlsx"),
  PANEL_Top2 = file.path(folder_etap02, "01_Tabele_Excel", "Tabela_09_PANEL_top2_peptydy_na_bialko_PEP.xlsx"),
  PANEL_Top3 = file.path(folder_etap02, "01_Tabele_Excel", "Tabela_07_PANEL_top3_peptydy_na_bialko_PEP.xlsx"),
  Status_PANEL = file.path(folder_etap02, "03_Kontrola", "Kontrola_PANEL_status_bialek_po_kryteriach_PEP.xlsx"),
  Podsumowanie_Top2_Top3 = file.path(folder_etap02, "04_CSV", "ETAP_02_porownanie_Top2_Top3.csv")
)

# ------------------------------------------------------------------------------
# 3. FUNKCJE
# ------------------------------------------------------------------------------

# Ujednolica wartości liczbowe wczytywane z arkuszy i plików CSV.
na_liczbe <- function(x) {
  if (is.numeric(x)) return(x)
  readr::parse_number(
    as.character(x),
    locale = readr::locale(decimal_mark = ",", grouping_mark = " ")
  )
}

czy_ragged <- function(peptyd) {
  stringr::str_detect(peptyd, wzorzec_ragged_N) |
    stringr::str_detect(peptyd, wzorzec_ragged_C)
}

# Przerywa wykonanie, jeżeli w tabeli brakuje kolumn potrzebnych w dalszej analizie.
sprawdz_kolumny <- function(df, wymagane, nazwa) {
  brak <- setdiff(wymagane, colnames(df))
  if (length(brak) > 0) {
    stop("W tabeli ", nazwa, " brakuje kolumn: ", paste(brak, collapse = ", "))
  }
}

zapisz_excel <- function(df, sciezka, arkusz = "Dane") {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, arkusz)
  openxlsx::writeData(wb, arkusz, df)
  openxlsx::freezePane(wb, arkusz, firstRow = TRUE)
  if (ncol(df) > 0) openxlsx::setColWidths(wb, arkusz, cols = seq_len(ncol(df)), widths = "auto")
  openxlsx::saveWorkbook(wb, sciezka, overwrite = TRUE)
}

zapisz_wiele <- function(lista, sciezka) {
  wb <- openxlsx::createWorkbook()
  for (nm in names(lista)) {
    ark <- substr(nm, 1, 31)
    openxlsx::addWorksheet(wb, ark)
    openxlsx::writeData(wb, ark, lista[[nm]])
    openxlsx::freezePane(wb, ark, firstRow = TRUE)
    if (ncol(lista[[nm]]) > 0) {
      openxlsx::setColWidths(wb, ark, cols = seq_len(ncol(lista[[nm]])), widths = "auto")
    }
  }
  openxlsx::saveWorkbook(wb, sciezka, overwrite = TRUE)
}

theme_prm <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2),
      axis.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(15, 25, 20, 20)
    )
}

zapisz_wykres <- function(p, nazwa, szer = 10, wys = 6) {
  ggplot2::ggsave(file.path(folder_wykresy, paste0(nazwa, ".png")), p, width = szer, height = wys, dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(folder_wykresy, paste0(nazwa, ".pdf")), p, width = szer, height = wys, bg = "white")
}

# Porządkuje listę Top 2/Top 3 i nadaje jednolity zestaw kolumn dla całej FASTA i panelu.
standaryzuj_liste <- function(df, zakres, typ_listy, maks_n) {
  wymagane <- c(
    "Bialko", "Nazwa_bialka", "Peptyd", "Wybrany_ladunek",
    "N_plikow_zmierzony", "RT_sredni", "RT_start_sredni", "RT_stop_sredni",
    "Najlepszy_PEP", "Max_intensywnosc", "Wybrane_mz",
    "Max_dobrych_fragmentow", "Max_suma_top_fragmentow",
    "Max_srednia_top_fragmentow", "Najmocniejszy_fragment",
    "Max_sredni_score_top3", "Max_min_score_top3", "Max_udzial_runow_top3",
    "Pozycja_w_bialku"
  )
  sprawdz_kolumny(df, wymagane, paste(zakres, typ_listy))

  kol_num <- setdiff(
    wymagane,
    c("Bialko", "Nazwa_bialka", "Peptyd")
  )

  df %>%
    select(-matches("QValue|Q[.]Value|Q_Value", ignore.case = TRUE)) %>%
    mutate(
      Bialko = as.character(Bialko),
      Nazwa_bialka = as.character(Nazwa_bialka),
      Peptyd = as.character(Peptyd)
    ) %>%
    mutate(across(all_of(kol_num), na_liczbe)) %>%
    mutate(
      Wybrany_ladunek = as.integer(Wybrany_ladunek),
      N_plikow_zmierzony = as.integer(N_plikow_zmierzony),
      Max_dobrych_fragmentow = as.integer(Max_dobrych_fragmentow),
      Pozycja_w_bialku = as.integer(Pozycja_w_bialku),
      Zakres_listy = zakres,
      Typ_listy = typ_listy,
      Maks_peptydow_na_bialko = maks_n,
      Ragged = czy_ragged(Peptyd),
      Zrodlo_RT = "obserwowany RT30",
      Zrodlo_fragmentow = "DIA-NN --export-quant",
      Okno_RT_empiryczne = RT_stop_sredni - RT_start_sredni,
      PEP_log10 = log10(pmax(Najlepszy_PEP, .Machine$double.xmin))
    ) %>%
    arrange(Bialko, Pozycja_w_bialku, Najlepszy_PEP, Peptyd)
}

# Wyznacza podstawowe liczebności gotowej listy peptydowo-prekursorowej.
podsumuj_liste <- function(df) {
  na_bialko <- df %>% count(Bialko, name = "Liczba_peptydow_na_bialko")
  tibble(
    Zakres_listy = unique(df$Zakres_listy),
    Typ_listy = unique(df$Typ_listy),
    Liczba_bialek = n_distinct(df$Bialko),
    Liczba_peptydow = nrow(df),
    Bialka_z_1_peptydem = sum(na_bialko$Liczba_peptydow_na_bialko == 1L),
    Bialka_z_2_peptydami = sum(na_bialko$Liczba_peptydow_na_bialko == 2L),
    Bialka_z_3_peptydami = sum(na_bialko$Liczba_peptydow_na_bialko == 3L),
    PEP_min = min(df$Najlepszy_PEP, na.rm = TRUE),
    PEP_mediana = median(df$Najlepszy_PEP, na.rm = TRUE),
    PEP_max = max(df$Najlepszy_PEP, na.rm = TRUE),
    Fragmenty_kwalifikujace_min = min(df$Max_dobrych_fragmentow, na.rm = TRUE),
    Fragmenty_kwalifikujace_mediana = median(df$Max_dobrych_fragmentow, na.rm = TRUE),
    Fragmenty_kwalifikujace_max = max(df$Max_dobrych_fragmentow, na.rm = TRUE),
    Score_top3_srednia = mean(df$Max_sredni_score_top3, na.rm = TRUE),
    Score_top3_mediana = median(df$Max_sredni_score_top3, na.rm = TRUE),
    Score_top3_min = min(df$Max_sredni_score_top3, na.rm = TRUE),
    RT_min = min(df$RT_sredni, na.rm = TRUE),
    RT_max = max(df$RT_sredni, na.rm = TRUE)
  )
}

# Zapisuje wynik pojedynczego testu poprawności.
kontrola <- function(nazwa, wynik, obserwowana = "", oczekiwana = "") {
  tibble(
    Kontrola = nazwa,
    Wynik = isTRUE(wynik),
    Wartosc_obserwowana = as.character(obserwowana),
    Wartosc_oczekiwana = as.character(oczekiwana)
  )
}

# ------------------------------------------------------------------------------
# 4. WCZYTANIE LIST I DYNAMICZNE OCZEKIWANE LICZEBNOŚCI Z ETAPU 02
# ------------------------------------------------------------------------------

# Liczebności kontrolne są pobierane z aktualnego etapu 02, a nie wpisywane ręcznie.

listy <- list(
  FASTA_Top2 = standaryzuj_liste(readxl::read_excel(sciezki$FASTA_Top2), "Cała FASTA", "Top 2", 2L),
  FASTA_Top3 = standaryzuj_liste(readxl::read_excel(sciezki$FASTA_Top3), "Cała FASTA", "Top 3", 3L),
  PANEL_Top2 = standaryzuj_liste(readxl::read_excel(sciezki$PANEL_Top2), "Panel 741 białek", "Top 2", 2L),
  PANEL_Top3 = standaryzuj_liste(readxl::read_excel(sciezki$PANEL_Top3), "Panel 741 białek", "Top 3", 3L)
)

podsumowanie_etap02 <- readr::read_csv2(
  sciezki$Podsumowanie_Top2_Top3,
  show_col_types = FALSE,
  locale = readr::locale(decimal_mark = ",", grouping_mark = " ", encoding = "UTF-8")
) %>%
  mutate(
    Limit_ID = as.integer(na_liczbe(Limit_ID)),
    Liczba_bialek = as.integer(na_liczbe(Liczba_bialek)),
    Liczba_peptydow = as.integer(na_liczbe(Liczba_peptydow))
  )

sprawdz_kolumny(
  podsumowanie_etap02,
  c("Zbior", "Limit_ID", "Liczba_bialek", "Liczba_peptydow"),
  "ETAP_02_porownanie_Top2_Top3.csv"
)

oczekiwane_dla <- function(zbior, limit_id, kolumna) {
  x <- podsumowanie_etap02 %>%
    filter(Zbior == zbior, Limit_ID == limit_id) %>%
    pull(all_of(kolumna))
  if (length(x) != 1L) stop("Nie udało się jednoznacznie odczytać wartości oczekiwanej z ETAPU 02.")
  x[[1]]
}

for (nm in names(listy)) {
  cat(
    nm, ": ", nrow(listy[[nm]]), " peptydów; ",
    n_distinct(listy[[nm]]$Bialko), " białek\n", sep = ""
  )
}

# ------------------------------------------------------------------------------
# 5. KONTROLE
# ------------------------------------------------------------------------------

# W tym miejscu sprawdzana jest zgodność limitu peptydów na białko, brak ragged oraz zgodność z wynikami etapu 02.

kontrole <- bind_rows(
  kontrola(
    "FASTA Top2 — liczba peptydów zgodna z ETAPEM 02",
    nrow(listy$FASTA_Top2) == oczekiwane_dla("Cała FASTA", 2L, "Liczba_peptydow"),
    nrow(listy$FASTA_Top2),
    oczekiwane_dla("Cała FASTA", 2L, "Liczba_peptydow")
  ),
  kontrola(
    "FASTA Top2 — liczba białek zgodna z ETAPEM 02",
    n_distinct(listy$FASTA_Top2$Bialko) == oczekiwane_dla("Cała FASTA", 2L, "Liczba_bialek"),
    n_distinct(listy$FASTA_Top2$Bialko),
    oczekiwane_dla("Cała FASTA", 2L, "Liczba_bialek")
  ),
  kontrola(
    "FASTA Top3 — liczba peptydów zgodna z ETAPEM 02",
    nrow(listy$FASTA_Top3) == oczekiwane_dla("Cała FASTA", 3L, "Liczba_peptydow"),
    nrow(listy$FASTA_Top3),
    oczekiwane_dla("Cała FASTA", 3L, "Liczba_peptydow")
  ),
  kontrola(
    "FASTA Top3 — liczba białek zgodna z ETAPEM 02",
    n_distinct(listy$FASTA_Top3$Bialko) == oczekiwane_dla("Cała FASTA", 3L, "Liczba_bialek"),
    n_distinct(listy$FASTA_Top3$Bialko),
    oczekiwane_dla("Cała FASTA", 3L, "Liczba_bialek")
  ),
  kontrola(
    "Panel Top2 — liczba peptydów zgodna z ETAPEM 02",
    nrow(listy$PANEL_Top2) == oczekiwane_dla("Panel 741 białek", 2L, "Liczba_peptydow"),
    nrow(listy$PANEL_Top2),
    oczekiwane_dla("Panel 741 białek", 2L, "Liczba_peptydow")
  ),
  kontrola(
    "Panel Top2 — liczba białek zgodna z ETAPEM 02",
    n_distinct(listy$PANEL_Top2$Bialko) == oczekiwane_dla("Panel 741 białek", 2L, "Liczba_bialek"),
    n_distinct(listy$PANEL_Top2$Bialko),
    oczekiwane_dla("Panel 741 białek", 2L, "Liczba_bialek")
  ),
  kontrola(
    "Panel Top3 — liczba peptydów zgodna z ETAPEM 02",
    nrow(listy$PANEL_Top3) == oczekiwane_dla("Panel 741 białek", 3L, "Liczba_peptydow"),
    nrow(listy$PANEL_Top3),
    oczekiwane_dla("Panel 741 białek", 3L, "Liczba_peptydow")
  ),
  kontrola(
    "Panel Top3 — liczba białek zgodna z ETAPEM 02",
    n_distinct(listy$PANEL_Top3$Bialko) == oczekiwane_dla("Panel 741 białek", 3L, "Liczba_bialek"),
    n_distinct(listy$PANEL_Top3$Bialko),
    oczekiwane_dla("Panel 741 białek", 3L, "Liczba_bialek")
  ),
  kontrola(
    "FASTA Top2 jest podzbiorem FASTA Top3",
    nrow(anti_join(listy$FASTA_Top2, listy$FASTA_Top3, by = c("Bialko", "Peptyd"))) == 0L,
    nrow(anti_join(listy$FASTA_Top2, listy$FASTA_Top3, by = c("Bialko", "Peptyd"))),
    0
  ),
  kontrola(
    "Panel Top2 jest podzbiorem Panel Top3",
    nrow(anti_join(listy$PANEL_Top2, listy$PANEL_Top3, by = c("Bialko", "Peptyd"))) == 0L,
    nrow(anti_join(listy$PANEL_Top2, listy$PANEL_Top3, by = c("Bialko", "Peptyd"))),
    0
  ),
  kontrola(
    "Wszystkie listy mają PEP <= 0,01",
    all(vapply(listy, function(x) all(!is.na(x$Najlepszy_PEP) & x$Najlepszy_PEP <= 0.01), logical(1))),
    "", "TRUE"
  ),
  kontrola(
    "Wszystkie listy mają >=3 kwalifikujące się fragmenty",
    all(vapply(listy, function(x) all(!is.na(x$Max_dobrych_fragmentow) & x$Max_dobrych_fragmentow >= 3L), logical(1))),
    "", "TRUE"
  ),
  kontrola(
    "Wszystkie listy mają obserwowany RT30",
    all(vapply(listy, function(x) all(is.finite(x$RT_sredni)), logical(1))),
    "", "TRUE"
  ),
  kontrola(
    "Wszystkie listy mają ładunki 2+ lub 3+",
    all(vapply(listy, function(x) all(x$Wybrany_ladunek %in% c(2L, 3L)), logical(1))),
    "", "TRUE"
  ),
  kontrola(
    "Wszystkie listy bez peptydów ragged",
    !any(vapply(listy, function(x) any(czy_ragged(x$Peptyd)), logical(1))),
    "", "TRUE"
  ),
  kontrola(
    "Wszystkie listy mają metrykę średniego Fr.Score top3",
    all(vapply(listy, function(x) all(is.finite(x$Max_sredni_score_top3)), logical(1))),
    "", "TRUE"
  ),
  kontrola(
    "Top2: maksymalnie dwa peptydy na białko",
    max(vapply(
      list(listy$FASTA_Top2, listy$PANEL_Top2),
      function(x) max(count(x, Bialko)$n),
      integer(1)
    )) <= 2L,
    max(vapply(
      list(listy$FASTA_Top2, listy$PANEL_Top2),
      function(x) max(count(x, Bialko)$n),
      integer(1)
    )),
    2
  ),
  kontrola(
    "Top3: maksymalnie trzy peptydy na białko",
    max(vapply(
      list(listy$FASTA_Top3, listy$PANEL_Top3),
      function(x) max(count(x, Bialko)$n),
      integer(1)
    )) <= 3L,
    max(vapply(
      list(listy$FASTA_Top3, listy$PANEL_Top3),
      function(x) max(count(x, Bialko)$n),
      integer(1)
    )),
    3
  )
)

zapisz_excel(kontrole, file.path(folder_kontrola, "KONTROLA_ETAP_03_FINAL.xlsx"), "Kontrole")
readr::write_excel_csv2(kontrole, file.path(folder_kontrola, "KONTROLA_ETAP_03_FINAL.csv"))

if (any(!kontrole$Wynik)) {
  print(kontrole %>% filter(!Wynik))
  stop("ETAP 03 zatrzymany: co najmniej jedna kontrola zakończyła się FALSE.")
}

# ------------------------------------------------------------------------------
# 6. PODSUMOWANIA I STATUS PANELU
# ------------------------------------------------------------------------------

podsumowanie_list <- bind_rows(lapply(listy, podsumuj_liste))

dystrybucja <- bind_rows(lapply(listy, function(df) {
  df %>%
    count(Bialko, name = "Liczba_peptydow_na_bialko") %>%
    count(Liczba_peptydow_na_bialko, name = "Liczba_bialek") %>%
    mutate(
      Zakres_listy = unique(df$Zakres_listy),
      Typ_listy = unique(df$Typ_listy)
    )
}))

status_panel <- readxl::read_excel(sciezki$Status_PANEL)
kol_status <- intersect(c("Status_koncowy", "Ostatni_etap"), colnames(status_panel))
if (length(kol_status) == 0L) {
  status_panel_podsumowanie <- tibble(
    Status = "Brak rozpoznanej kolumny statusu",
    Liczba_bialek = nrow(status_panel)
  )
} else {
  status_panel_podsumowanie <- status_panel %>%
    group_by(across(all_of(kol_status[[1]]))) %>%
    summarise(Liczba_bialek = n(), .groups = "drop") %>%
    rename(Status = all_of(kol_status[[1]])) %>%
    mutate(Procent_panelu = round(100 * Liczba_bialek / sum(Liczba_bialek), 2)) %>%
    arrange(desc(Liczba_bialek))
}

opis_metody <- tibble::tribble(
  ~Element, ~Opis,
  "Źródło", "Finalny ETAP 02 oparty na głównych raportach DIA-NN z --export-quant.",
  "PEP", "PEP <= 0,01; Q.Value nie uczestniczy w selekcji ani rankingu.",
  "Fragmenty", "Minimum 3 kwalifikujące się fragmenty. Fr.Score nie ma arbitralnego progu i służy do rankingu top 3 fragmentów.",
  "RT", "Bezpośrednio obserwowany RT w gradiencie 30-minutowym.",
  "Listy", "Top 2 i Top 3 dla całej FASTA oraz panelu 741 białek.",
  "Charakter wyniku", "Listy peptydowo-prekursorowe do dalszego planowania PRM i oceny w Skyline."
)

zapisz_wiele(
  list(
    Metoda = opis_metody,
    Podsumowanie_list = podsumowanie_list,
    Dystrybucja = dystrybucja,
    Status_PANEL = status_panel_podsumowanie,
    Kontrole = kontrole
  ),
  file.path(folder_tabele, "Tabela_00_Podsumowanie_ETAP_03_EXPORT_QUANT.xlsx")
)

# ------------------------------------------------------------------------------
# 7. EKSPORT LIST
# ------------------------------------------------------------------------------

# Do dalszych etapów zapisywane są oddzielne listy Top 2 i Top 3 dla całej FASTA oraz panelu.

nazwy_xlsx <- c(
  FASTA_Top2 = "Tabela_01A_Lista_peptydowo_prekursorowa_CALA_FASTA_top2_PEP_RT30_PO_RAGGED.xlsx",
  FASTA_Top3 = "Tabela_01B_Lista_peptydowo_prekursorowa_CALA_FASTA_top3_PEP_RT30_PO_RAGGED.xlsx",
  PANEL_Top2 = "Tabela_02A_Lista_peptydowo_prekursorowa_PANEL_741_top2_PEP_RT30_PO_RAGGED.xlsx",
  PANEL_Top3 = "Tabela_02B_Lista_peptydowo_prekursorowa_PANEL_741_top3_PEP_RT30_PO_RAGGED.xlsx"
)

for (nm in names(listy)) {
  zapisz_excel(listy[[nm]], file.path(folder_tabele, nazwy_xlsx[[nm]]), "Lista")
}

# Aliasy zgodnościowe starszego ETAPU 03 — odpowiadają Top3.
zapisz_excel(
  listy$FASTA_Top3,
  file.path(folder_tabele, "Tabela_01_Lista_pomiarowa_CALA_FASTA_top3_PEP_RT30.xlsx"),
  "Lista"
)
zapisz_excel(
  listy$PANEL_Top3,
  file.path(folder_tabele, "Tabela_02_Lista_pomiarowa_PANEL_741_top3_PEP_RT30.xlsx"),
  "Lista"
)

nazwy_csv <- c(
  FASTA_Top2 = "CSV_01A_CALA_FASTA_lista_peptydowo_prekursorowa_top2_PO_RAGGED.csv",
  FASTA_Top3 = "CSV_01B_CALA_FASTA_lista_peptydowo_prekursorowa_top3_PO_RAGGED.csv",
  PANEL_Top2 = "CSV_02A_PANEL_741_lista_peptydowo_prekursorowa_top2_PO_RAGGED.csv",
  PANEL_Top3 = "CSV_02B_PANEL_741_lista_peptydowo_prekursorowa_top3_PO_RAGGED.csv"
)
for (nm in names(listy)) {
  readr::write_excel_csv2(listy[[nm]], file.path(folder_csv, nazwy_csv[[nm]]))
}

# Aliasy zgodnościowe odpowiadające wcześniejszemu ETAPOWI 03 (Top 3).
readr::write_excel_csv2(
  listy$FASTA_Top3,
  file.path(folder_csv, "CSV_01_CALA_FASTA_lista_peptydowo_prekursorowa_PO_RAGGED.csv")
)
readr::write_excel_csv2(
  listy$PANEL_Top3,
  file.path(folder_csv, "CSV_02_PANEL_741_lista_peptydowo_prekursorowa_PO_RAGGED.csv")
)

# ------------------------------------------------------------------------------
# 8. WYKRESY
# ------------------------------------------------------------------------------

p1 <- podsumowanie_list %>%
  select(Zakres_listy, Typ_listy, Liczba_bialek, Liczba_peptydow) %>%
  pivot_longer(c(Liczba_bialek, Liczba_peptydow), names_to = "Miara", values_to = "Liczba") %>%
  mutate(Miara = recode(Miara, Liczba_bialek = "Białka", Liczba_peptydow = "Peptydy")) %>%
  ggplot(aes(x = Typ_listy, y = Liczba, fill = Miara)) +
  geom_col(position = "dodge") +
  facet_wrap(~Zakres_listy, scales = "free_y") +
  labs(
    title = "Listy peptydowo-prekursorowe Top 2 i Top 3",
    subtitle = "PEP <= 0,01; >=3 kwalifikujące się fragmenty; obserwowany RT30",
    x = NULL, y = "Liczba", fill = "Miara"
  ) +
  theme_prm()
zapisz_wykres(p1, "Wykres_01_Liczba_bialek_i_peptydow_Top2_Top3", 11, 7)

p2 <- bind_rows(listy) %>%
  ggplot(aes(x = Max_dobrych_fragmentow)) +
  geom_histogram(binwidth = 1, boundary = 0) +
  facet_grid(Zakres_listy ~ Typ_listy, scales = "free_y") +
  labs(
    title = "Liczba kwalifikujących się fragmentów",
    subtitle = "Fragmenty ocenione na podstawie danych DIA-NN --export-quant",
    x = "Liczba kwalifikujących się fragmentów",
    y = "Liczba peptydów"
  ) +
  theme_prm()
zapisz_wykres(p2, "Wykres_02_Rozklad_kwalifikujacych_fragmentow", 11, 7)

p3 <- bind_rows(listy) %>%
  ggplot(aes(x = Max_sredni_score_top3)) +
  geom_histogram(bins = 40) +
  facet_grid(Zakres_listy ~ Typ_listy, scales = "free_y") +
  labs(
    title = "Rozkład średniego Fr.Score dla trzech najlepszych fragmentów",
    subtitle = "Brak arbitralnego progu Fr.Score; metryka służy do rankingu jakości fragmentów",
    x = "Średni Fr.Score top 3",
    y = "Liczba peptydów"
  ) +
  theme_prm()
zapisz_wykres(p3, "Wykres_03_Rozklad_FrScore_top3", 11, 7)

# ------------------------------------------------------------------------------
# 9. INDEKS I RAPORT
# ------------------------------------------------------------------------------

indeks <- tibble::tribble(
  ~Folder, ~Plik, ~Opis,
  "01_Tabele_Excel", "Tabela_00_Podsumowanie_ETAP_03_EXPORT_QUANT.xlsx", "Podsumowanie, status panelu i kontrole.",
  "01_Tabele_Excel", nazwy_xlsx[["FASTA_Top2"]], "Cała FASTA — Top 2.",
  "01_Tabele_Excel", nazwy_xlsx[["FASTA_Top3"]], "Cała FASTA — Top 3.",
  "01_Tabele_Excel", nazwy_xlsx[["PANEL_Top2"]], "Panel 741 — Top 2.",
  "01_Tabele_Excel", nazwy_xlsx[["PANEL_Top3"]], "Panel 741 — Top 3.",
  "02_CSV_do_dalszego_przygotowania", nazwy_csv[["FASTA_Top2"]], "Wejście Top2 FASTA dla ETAPU 04.",
  "02_CSV_do_dalszego_przygotowania", nazwy_csv[["FASTA_Top3"]], "Wejście Top3 FASTA dla ETAPU 04.",
  "02_CSV_do_dalszego_przygotowania", nazwy_csv[["PANEL_Top2"]], "Wejście Top2 panel dla ETAPÓW 04/05.",
  "02_CSV_do_dalszego_przygotowania", nazwy_csv[["PANEL_Top3"]], "Wejście Top3 panel dla ETAPÓW 04/05.",
  "04_Kontrola", "KONTROLA_ETAP_03_FINAL.xlsx/csv", "Kontrole zgodności z finalnym ETAPEM 02 i warunków list."
)
zapisz_excel(indeks, file.path(folder_wyniki, "00_INDEKS_PLIKOW_ETAP_03_EXPORT_QUANT.xlsx"), "Indeks")

raport <- capture.output({
  cat("ETAP 03 — LISTY PEPTYDOWO-PREKURSOROWE — EXPORT-QUANT\n")
  cat("Data:", format(Sys.time()), "\n")
  cat("Źródło ETAP 02:", folder_etap02, "\n\n")
  cat("PODSUMOWANIE\n")
  print(podsumowanie_list)
  cat("\nKONTROLE\n")
  print(kontrole)
  cat("\nUwagi:\n")
  cat("- PEP <= 0,01 jest kryterium jakości identyfikacji.\n")
  cat("- Minimum 3 kwalifikujące się fragmenty pochodzi z danych EXPORT-QUANT.\n")
  cat("- Fr.Score nie ma arbitralnego progu; zachowano go jako metrykę rankingu.\n")
  cat("- RT pochodzi bezpośrednio z gradientu 30-minutowego.\n")
})
writeLines(raport, file.path(folder_kontrola, "RAPORT_ETAP_03_FINAL.txt"), useBytes = TRUE)

cat("\nGOTOWE — ETAP 03 EXPORT-QUANT\n")
cat("Folder wynikowy:\n", normalizePath(folder_wyniki, winslash = "/"), "\n")
print(podsumowanie_list)

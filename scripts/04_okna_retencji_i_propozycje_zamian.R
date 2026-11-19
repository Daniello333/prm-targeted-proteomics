# ==============================================================================
# 04. ANALIZA ZATŁOCZENIA CZASU RETENCJI I PROPOZYCJE ZAMIAN
#
# Praca magisterska:
# „Planowanie celowanych eksperymentów proteomicznych z wykorzystaniem środowiska R”
# Autor: Daniel Drabik
#
# Skrypt ocenia rozmieszczenie peptydów z list Top 2 i Top 3 wzdłuż gradientu.
# Zatłoczenie jest liczone dla stałych okien 1-minutowych, ruchomego otoczenia
# RT ±1 min oraz na podstawie nakładania empirycznych przedziałów RT.Start–RT.Stop.
#
# Dla list panelowych sprawdzane są także możliwe zamiany peptydów w obrębie
# tego samego białka. Zamiana jest traktowana wyłącznie jako propozycja do
# dalszej oceny; skrypt nie uznaje jej automatycznie za finalny wybór do PRM.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. PAKIETY
# ------------------------------------------------------------------------------

pakiety <- c("tidyverse", "readxl", "openxlsx", "ggplot2")
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

folder_etap03 <- NULL
folder_etap02 <- NULL

szerokosc_okna_stalego_min <- 1
polowa_okna_ruchomego_min <- 1

prog_PEP <- 0.01
min_dobrych_fragmentow <- 3L

# Poziom TOLERANCJA — parametry „podobnej jakości”.
maksymalny_mnoznik_PEP <- 2
maksymalny_bezwzgledny_wzrost_PEP <- 0.0005
minimalny_udzial_intensywnosci <- 0.90

# Minimalna lokalna poprawa zatłoczenia wymagana do rozważenia ruchu.
minimalna_poprawa_pm1 <- 5L
minimalna_poprawa_okno_1min <- 1L
minimalna_poprawa_empiryczna <- 1L

chron_bialka_z_jednym_peptydem <- TRUE
maksymalna_liczba_zamian_na_bialko <- 1L

# Liczebności kontrolne nie są wpisane na sztywno. Zostaną odczytane
# bezpośrednio z plików podsumowujących finalnego ETAPU 02.

znacznik <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Ustala katalog skryptu i lokalizację nowych wyników.
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
folder_wyniki <- file.path(
  katalog_skryptu,
  paste0(
    "WYNIKI_ETAP_04_OKNA_RETENCJI_I_PROPOZYCJE_ZAMIAN_EXPORT_QUANT_HUMAN_ONLY_",
    znacznik
  )
)
folder_tabele <- file.path(folder_wyniki, "01_Tabele_Excel")
folder_wykresy <- file.path(folder_wyniki, "02_Wykresy")
folder_kontrola <- file.path(folder_wyniki, "03_Kontrola")
folder_csv <- file.path(folder_wyniki, "04_CSV")
for (f in c(folder_wyniki, folder_tabele, folder_wykresy, folder_kontrola, folder_csv)) {
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
# 2. WYSZUKIWANIE FOLDERÓW LUB ZIP-ÓW WEJŚCIOWYCH
# ------------------------------------------------------------------------------

# Sprawdza, czy katalog zawiera komplet wymaganych plików z wcześniejszego etapu.
ma_pliki <- function(folder, wymagane) {
  dir.exists(folder) && all(file.exists(file.path(folder, wymagane)))
}

# Wybiera najnowszą zgodną wersję wyników.
wybierz_najnowszy <- function(sciezki) {
  sciezki <- unique(sciezki[file.exists(sciezki) | dir.exists(sciezki)])
  if (length(sciezki) == 0) return(NULL)
  sciezki[which.max(file.info(sciezki)$mtime)]
}

lista_folderow <- function() {
  unique(unlist(lapply(korzenie_szukania, function(root) {
    c(root, list.dirs(root, recursive = TRUE, full.names = TRUE))
  })))
}

lista_zipow <- function() {
  unique(unlist(lapply(korzenie_szukania, function(root) {
    list.files(
      root,
      pattern = "[.]zip$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )
  })))
}

zip_zawiera <- function(zip_path, wymagane) {
  listing <- tryCatch(utils::unzip(zip_path, list = TRUE), error = function(e) NULL)
  if (is.null(listing)) return(FALSE)
  nazwy <- gsub("\\\\", "/", listing$Name)
  all(vapply(wymagane, function(rel) {
    rel2 <- gsub("\\\\", "/", rel)
    any(endsWith(nazwy, rel2))
  }, logical(1)))
}

wypakuj_do_temp <- function(zip_path, wymagane, etykieta) {
  katalog_temp <- tempfile(pattern = paste0("ETAP04_", etykieta, "_"))
  dir.create(katalog_temp, recursive = TRUE, showWarnings = FALSE)
  utils::unzip(zip_path, exdir = katalog_temp)

  kandydaci <- c(katalog_temp, list.dirs(katalog_temp, recursive = TRUE, full.names = TRUE))
  kandydaci <- kandydaci[vapply(kandydaci, ma_pliki, logical(1), wymagane = wymagane)]
  wybierz_najnowszy(kandydaci)
}

# Szuka wyników wcześniejszego etapu w katalogach lub archiwach ZIP.
znajdz_lub_wypakuj <- function(wymagane, preferowany_wzorzec, wzorzec_zip, etykieta) {
  foldery <- lista_folderow()
  kandydaci <- foldery[vapply(foldery, ma_pliki, logical(1), wymagane = wymagane)]

  preferowane <- kandydaci[stringr::str_detect(
    basename(kandydaci),
    stringr::regex(preferowany_wzorzec, ignore_case = TRUE)
  )]
  if (length(preferowane) > 0) return(wybierz_najnowszy(preferowane))
  if (length(kandydaci) > 0) return(wybierz_najnowszy(kandydaci))

  zipy <- lista_zipow()
  zipy <- zipy[stringr::str_detect(
    basename(zipy),
    stringr::regex(wzorzec_zip, ignore_case = TRUE)
  )]
  zipy <- zipy[vapply(zipy, zip_zawiera, logical(1), wymagane = wymagane)]
  if (length(zipy) == 0) return(NULL)

  zip_najnowszy <- wybierz_najnowszy(zipy)
  wypakuj_do_temp(zip_najnowszy, wymagane, etykieta)
}

wymagane_etap03 <- c(
  file.path(
    "02_CSV_do_dalszego_przygotowania",
    "CSV_01A_CALA_FASTA_lista_peptydowo_prekursorowa_top2_PO_RAGGED.csv"
  ),
  file.path(
    "02_CSV_do_dalszego_przygotowania",
    "CSV_01B_CALA_FASTA_lista_peptydowo_prekursorowa_top3_PO_RAGGED.csv"
  ),
  file.path(
    "02_CSV_do_dalszego_przygotowania",
    "CSV_02A_PANEL_741_lista_peptydowo_prekursorowa_top2_PO_RAGGED.csv"
  ),
  file.path(
    "02_CSV_do_dalszego_przygotowania",
    "CSV_02B_PANEL_741_lista_peptydowo_prekursorowa_top3_PO_RAGGED.csv"
  )
)

wymagane_etap02 <- c(
  file.path(
    "01_Tabele_Excel",
    "Tabela_05_PANEL_kandydaci_po_PEP_0_01_i_min_3_fragmentach.xlsx"
  ),
  file.path(
    "04_CSV",
    "ETAP_02_porownanie_Top2_Top3.csv"
  ),
  file.path(
    "04_CSV",
    "ETAP_02_redukcja_PANEL_741.csv"
  )
)

if (is.null(folder_etap03)) {
  folder_etap03 <- znajdz_lub_wypakuj(
    wymagane_etap03,
    "^WYNIKI_ETAP_03_.*HUMAN_ONLY_",
    "^WYNIKI_ETAP_03_.*HUMAN_ONLY_.*[.]zip$",
    "ETAP03"
  )
}
if (is.null(folder_etap02)) {
  folder_etap02 <- znajdz_lub_wypakuj(
    wymagane_etap02,
    "^WYNIKI_ETAP_02_.*EXPORT_QUANT_HUMAN_ONLY_",
    "^WYNIKI_ETAP_02_.*EXPORT_QUANT_HUMAN_ONLY_.*[.]zip$",
    "ETAP02"
  )
}

if (is.null(folder_etap03)) {
  stop("Nie znaleziono finalnego ETAPU 03 HUMAN-ONLY. Ustaw folder_etap03 ręcznie.")
}
if (is.null(folder_etap02)) {
  stop("Nie znaleziono finalnego ETAPU 02 HUMAN-ONLY. Ustaw folder_etap02 ręcznie.")
}

folder_etap03 <- normalizePath(folder_etap03, winslash = "/", mustWork = TRUE)
folder_etap02 <- normalizePath(folder_etap02, winslash = "/", mustWork = TRUE)

cat("ETAP 03:\n", folder_etap03, "\n\n")
cat("ETAP 02:\n", folder_etap02, "\n\n")
cat("Wyniki ETAPU 04:\n", folder_wyniki, "\n\n")


# ------------------------------------------------------------------------------
# 3. FUNKCJE POMOCNICZE
# ------------------------------------------------------------------------------

# Ujednolica zapis wartości liczbowych po wczytaniu danych.
na_liczbe <- function(x) {
  if (is.numeric(x)) return(x)
  readr::parse_number(
    as.character(x),
    locale = readr::locale(decimal_mark = ",", grouping_mark = " ")
  )
}

czy_ragged <- function(peptyd) {
  stringr::str_detect(peptyd, "^(KK|KR|RK|RR)") |
    stringr::str_detect(peptyd, "(KK|KR|RK|RR)$")
}

sprawdz_kolumny <- function(df, wymagane, nazwa) {
  brak <- setdiff(wymagane, colnames(df))
  if (length(brak) > 0) {
    stop("W tabeli ", nazwa, " brakuje kolumn: ", paste(brak, collapse = ", "))
  }
}

# Porządkuje bazową listę Top 2/Top 3 i sprawdza kolumny potrzebne do analizy czasu retencji.
przygotuj_liste <- function(df, zakres, typ_listy) {
  wymagane <- c(
    "Bialko", "Nazwa_bialka", "Peptyd", "Wybrany_ladunek",
    "N_plikow_zmierzony", "RT_sredni", "RT_start_sredni", "RT_stop_sredni",
    "Najlepszy_PEP", "Max_intensywnosc", "Wybrane_mz",
    "Max_dobrych_fragmentow", "Max_srednia_top_fragmentow",
    "Max_suma_top_fragmentow", "Najmocniejszy_fragment",
    "Max_sredni_score_top3", "Max_min_score_top3", "Max_udzial_runow_top3",
    "Pozycja_w_bialku"
  )
  sprawdz_kolumny(df, wymagane, paste(zakres, typ_listy))

  kol_num <- setdiff(wymagane, c("Bialko", "Nazwa_bialka", "Peptyd"))

  df %>%
    mutate(
      Bialko = as.character(Bialko),
      Nazwa_bialka = as.character(Nazwa_bialka),
      Peptyd = as.character(Peptyd)
    ) %>%
    mutate(across(all_of(kol_num), na_liczbe)) %>%
    transmute(
      Bialko,
      Nazwa_bialka,
      Peptyd,
      Wybrany_ladunek = as.integer(Wybrany_ladunek),
      N_plikow_zmierzony = as.integer(N_plikow_zmierzony),
      RT_sredni,
      RT_start_sredni,
      RT_stop_sredni,
      Najlepszy_PEP,
      Max_intensywnosc,
      Wybrane_mz,
      Max_dobrych_fragmentow = as.integer(Max_dobrych_fragmentow),
      Max_srednia_top_fragmentow,
      Max_suma_top_fragmentow,
      Najmocniejszy_fragment,
      Max_sredni_score_top3,
      Max_min_score_top3,
      Max_udzial_runow_top3,
      Pozycja_w_bialku = as.integer(Pozycja_w_bialku),
      Zakres_analizy = zakres,
      Typ_listy = typ_listy,
      Klucz = paste(Bialko, Peptyd, sep = "||")
    ) %>%
    distinct(Bialko, Peptyd, .keep_all = TRUE)
}

# Przygotowuje pełną pulę panelowych kandydatów, z której mogą być proponowane zamiany.
przygotuj_kandydatow <- function(df) {
  wymagane <- c(
    "Bialko", "Nazwa_bialka", "Peptyd", "Wybrany_ladunek",
    "N_plikow_zmierzony", "RT_sredni", "RT_start_sredni", "RT_stop_sredni",
    "Najlepszy_PEP", "Max_intensywnosc", "Wybrane_mz",
    "Max_dobrych_fragmentow", "Max_srednia_top_fragmentow",
    "Max_suma_top_fragmentow", "Najmocniejszy_fragment",
    "Max_sredni_score_top3", "Max_min_score_top3", "Max_udzial_runow_top3"
  )
  sprawdz_kolumny(df, wymagane, "pełna pula kandydatów panelu")

  kol_num <- setdiff(wymagane, c("Bialko", "Nazwa_bialka", "Peptyd"))

  df %>%
    mutate(
      Bialko = as.character(Bialko),
      Nazwa_bialka = as.character(Nazwa_bialka),
      Peptyd = as.character(Peptyd)
    ) %>%
    mutate(across(all_of(kol_num), na_liczbe)) %>%
    transmute(
      Bialko,
      Nazwa_bialka,
      Peptyd,
      Wybrany_ladunek = as.integer(Wybrany_ladunek),
      N_plikow_zmierzony = as.integer(N_plikow_zmierzony),
      RT_sredni,
      RT_start_sredni,
      RT_stop_sredni,
      Najlepszy_PEP,
      Max_intensywnosc,
      Wybrane_mz,
      Max_dobrych_fragmentow = as.integer(Max_dobrych_fragmentow),
      Max_srednia_top_fragmentow,
      Max_suma_top_fragmentow,
      Najmocniejszy_fragment,
      Max_sredni_score_top3,
      Max_min_score_top3,
      Max_udzial_runow_top3,
      Klucz = paste(Bialko, Peptyd, sep = "||")
    ) %>%
    distinct(Bialko, Peptyd, .keep_all = TRUE) %>%
    mutate(
      .rank_N_plikow = replace_na(N_plikow_zmierzony, 0L),
      .rank_PEP = replace_na(Najlepszy_PEP, Inf),
      .rank_score = replace_na(Max_sredni_score_top3, -Inf),
      .rank_suma = replace_na(Max_suma_top_fragmentow, 0),
      .rank_srednia = replace_na(Max_srednia_top_fragmentow, 0),
      .rank_int = replace_na(Max_intensywnosc, 0),
      .rank_fragmenty = replace_na(Max_dobrych_fragmentow, 0L)
    ) %>%
    group_by(Bialko) %>%
    arrange(
      desc(.rank_N_plikow),
      .rank_PEP,
      desc(.rank_score),
      desc(.rank_suma),
      desc(.rank_srednia),
      desc(.rank_int),
      desc(.rank_fragmenty),
      Peptyd,
      .by_group = TRUE
    ) %>%
    mutate(
      Pelny_ranking_w_bialku = row_number(),
      Liczba_kandydatow_w_bialku = n()
    ) %>%
    ungroup() %>%
    select(-starts_with(".rank_"))
}

czytaj_csv2 <- function(path) {
  readr::read_csv2(path, show_col_types = FALSE, progress = FALSE)
}

zapisz_excel <- function(df, path, arkusz = "Dane") {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, arkusz)
  openxlsx::writeData(wb, arkusz, df)
  openxlsx::freezePane(wb, arkusz, firstRow = TRUE)
  if (ncol(df) > 0) {
    openxlsx::setColWidths(wb, arkusz, cols = seq_len(ncol(df)), widths = "auto")
  }
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
}

zapisz_wiele_arkuszy <- function(lista_df, path) {
  wb <- openxlsx::createWorkbook()
  for (nm in names(lista_df)) {
    arkusz <- substr(nm, 1, 31)
    openxlsx::addWorksheet(wb, arkusz)
    openxlsx::writeData(wb, arkusz, lista_df[[nm]])
    openxlsx::freezePane(wb, arkusz, firstRow = TRUE)
    if (ncol(lista_df[[nm]]) > 0) {
      openxlsx::setColWidths(
        wb, arkusz, cols = seq_len(ncol(lista_df[[nm]])), widths = "auto"
      )
    }
  }
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
}

zapisz_csv2 <- function(df, nazwa) {
  readr::write_excel_csv2(df, file.path(folder_csv, nazwa))
}

# Dla każdego RT liczy cele znajdujące się w zadanym ruchomym otoczeniu czasowym.
licz_punkty_wokol <- function(punkty, wszystkie_rt, half_width = 1) {
  if (length(wszystkie_rt) == 0) return(integer(length(punkty)))
  posortowane <- sort(wszystkie_rt)
  eps <- 1e-10
  gorne <- findInterval(punkty + half_width + eps, posortowane)
  dolne <- findInterval(punkty - half_width - eps, posortowane)
  pmax(gorne - dolne, 0L)
}

# Dla każdego peptydu liczy liczbę nakładających się przedziałów RT.Start–RT.Stop.
licz_nakladanie_empiryczne <- function(starty, stopy) {
  posortowane_starty <- sort(starty)
  posortowane_stopy <- sort(stopy)
  eps <- 1e-10
  findInterval(stopy + eps, posortowane_starty) -
    findInterval(starty - eps, posortowane_stopy)
}

# Przypisuje peptydy do stałych okien 1-minutowych i oblicza ich liczebność.
analiza_okien_1min <- function(df) {
  szer <- szerokosc_okna_stalego_min
  df2 <- df %>%
    mutate(
      Okno_start = floor(RT_sredni / szer) * szer,
      Okno_koniec = Okno_start + szer
    )

  min_start <- floor(min(df2$RT_sredni) / szer) * szer
  max_start <- floor(max(df2$RT_sredni) / szer) * szer
  wszystkie <- tibble(Okno_start = seq(min_start, max_start, by = szer)) %>%
    mutate(Okno_koniec = Okno_start + szer)

  df2 %>%
    count(Okno_start, Okno_koniec, name = "Liczba_peptydow") %>%
    right_join(wszystkie, by = c("Okno_start", "Okno_koniec")) %>%
    mutate(
      Liczba_peptydow = replace_na(Liczba_peptydow, 0L),
      Zakres_analizy = unique(df$Zakres_analizy),
      Typ_listy = unique(df$Typ_listy),
      Okno_RT = paste0(Okno_start, "–", Okno_koniec, " min")
    ) %>%
    arrange(Okno_start) %>%
    select(Zakres_analizy, Typ_listy, everything())
}

# Oblicza zatłoczenie w ruchomym otoczeniu RT ±1 min.
analiza_pm1 <- function(df) {
  df %>%
    mutate(
      Liczba_wspolwystepujacych_pm1 = licz_punkty_wokol(
        RT_sredni, RT_sredni, polowa_okna_ruchomego_min
      ),
      Okno_pm1_start = pmax(RT_sredni - polowa_okna_ruchomego_min, 0),
      Okno_pm1_koniec = RT_sredni + polowa_okna_ruchomego_min
    )
}

# Oblicza zatłoczenie wynikające z nakładania empirycznych przedziałów retencji.
analiza_empiryczna <- function(df) {
  df %>%
    mutate(
      Liczba_nakladajacych_empirycznych = licz_nakladanie_empiryczne(
        RT_start_sredni, RT_stop_sredni
      ),
      Szerokosc_empirycznego_okna = RT_stop_sredni - RT_start_sredni
    )
}

podsumuj_zatloczenie <- function(df, wariant = "Lista bazowa") {
  okna <- analiza_okien_1min(df)
  pm <- analiza_pm1(df)
  emp <- analiza_empiryczna(df)

  tibble(
    Zakres_analizy = unique(df$Zakres_analizy),
    Typ_listy = unique(df$Typ_listy),
    Wariant = wariant,
    Liczba_bialek = n_distinct(df$Bialko),
    Liczba_peptydow = nrow(df),
    Maks_okno_1min = max(okna$Liczba_peptydow),
    Maks_pm1 = max(pm$Liczba_wspolwystepujacych_pm1),
    Maks_empiryczne = max(emp$Liczba_nakladajacych_empirycznych),
    Suma_kwadratow_okien_1min = sum(okna$Liczba_peptydow^2),
    Suma_pm1 = sum(pm$Liczba_wspolwystepujacych_pm1),
    Suma_empiryczne = sum(emp$Liczba_nakladajacych_empirycznych)
  )
}

# Dołącza do listy trzy lokalne miary zatłoczenia używane przy ocenie zamian.
dodaj_lokalne_zatloczenie <- function(df) {
  okna <- analiza_okien_1min(df) %>%
    select(Okno_start, Liczba_peptydow)

  analiza_pm1(df) %>%
    analiza_empiryczna() %>%
    mutate(Okno_start = floor(RT_sredni / szerokosc_okna_stalego_min) *
             szerokosc_okna_stalego_min) %>%
    left_join(okna, by = "Okno_start") %>%
    rename(Liczba_w_oknie_1min = Liczba_peptydow)
}

pusta_tabela_ruchow <- function() {
  tibble(
    Poziom_jakosci = character(),
    Bialko = character(),
    Nazwa_bialka = character(),
    Stary_peptyd = character(),
    Nowy_peptyd = character(),
    Stary_RT = numeric(),
    Nowy_RT = numeric(),
    Stary_PEP = numeric(),
    Nowy_PEP = numeric(),
    Mnoznik_PEP = numeric(),
    Bezwzgledny_wzrost_PEP = numeric(),
    Stare_N_plikow = integer(),
    Nowe_N_plikow = integer(),
    Stare_fragmenty = integer(),
    Nowe_fragmenty = integer(),
    Stary_sredni_score_top3 = numeric(),
    Nowy_sredni_score_top3 = numeric(),
    Roznica_sredniego_score_top3 = numeric(),
    Stary_min_score_top3 = numeric(),
    Nowy_min_score_top3 = numeric(),
    Roznica_min_score_top3 = numeric(),
    Stary_udzial_runow_top3 = numeric(),
    Nowy_udzial_runow_top3 = numeric(),
    Roznica_udzialu_runow_top3 = numeric(),
    Udzial_sumy_fragmentow = numeric(),
    Udzial_sredniej_fragmentow = numeric(),
    Udzial_max_intensywnosci = numeric(),
    Stare_okno_1min = integer(),
    Nowe_okno_1min = integer(),
    Poprawa_okno_1min = integer(),
    Stare_pm1 = integer(),
    Nowe_pm1 = integer(),
    Poprawa_pm1 = integer(),
    Stare_empiryczne = integer(),
    Nowe_empiryczne = integer(),
    Poprawa_empiryczna = integer(),
    Nowy_pelny_ranking = integer()
  )
}

# Tworzy dopuszczalne propozycje zamiany peptydu na alternatywnego kandydata tego samego białka.
enumeruj_ruchy <- function(
    aktualna,
    kandydaci,
    poziom = c("STRICT", "TOLERANCJA"),
    uzyte_bialka = character(0),
    wycofane_klucze = character(0)
) {
  poziom <- match.arg(poziom)
  lokalne <- dodaj_lokalne_zatloczenie(aktualna)
  liczby_bialek <- aktualna %>% count(Bialko, name = "Liczba_peptydow_w_bialku")
  lokalne <- lokalne %>% left_join(liczby_bialek, by = "Bialko")

  wybrane_klucze <- aktualna$Klucz
  rt_biezace <- aktualna$RT_sredni
  starty_biezace <- aktualna$RT_start_sredni
  stopy_biezace <- aktualna$RT_stop_sredni
  okna_biezace <- floor(rt_biezace / szerokosc_okna_stalego_min) *
    szerokosc_okna_stalego_min
  tab_okien <- table(okna_biezace)

  wyniki <- list()

  for (i in seq_len(nrow(lokalne))) {
    stary <- lokalne[i, ]
    if (stary$Bialko[[1]] %in% uzyte_bialka) next
    if (
      chron_bialka_z_jednym_peptydem &&
      stary$Liczba_peptydow_w_bialku[[1]] == 1L
    ) next

    alternatywy <- kandydaci %>%
      filter(Bialko == stary$Bialko[[1]]) %>%
      filter(!Klucz %in% wybrane_klucze, !Klucz %in% wycofane_klucze)
    if (nrow(alternatywy) == 0) next

    if (poziom == "STRICT") {
      alternatywy <- alternatywy %>%
        filter(
          N_plikow_zmierzony >= stary$N_plikow_zmierzony[[1]],
          Najlepszy_PEP <= stary$Najlepszy_PEP[[1]],
          Max_dobrych_fragmentow >= stary$Max_dobrych_fragmentow[[1]],
          Max_sredni_score_top3 >= stary$Max_sredni_score_top3[[1]],
          Max_min_score_top3 >= stary$Max_min_score_top3[[1]],
          Max_udzial_runow_top3 >= stary$Max_udzial_runow_top3[[1]],
          Max_suma_top_fragmentow >= stary$Max_suma_top_fragmentow[[1]],
          Max_srednia_top_fragmentow >= stary$Max_srednia_top_fragmentow[[1]],
          Max_intensywnosc >= stary$Max_intensywnosc[[1]]
        )
    } else {
      limit_PEP <- min(
        stary$Najlepszy_PEP[[1]] * maksymalny_mnoznik_PEP,
        stary$Najlepszy_PEP[[1]] + maksymalny_bezwzgledny_wzrost_PEP
      )
      alternatywy <- alternatywy %>%
        filter(
          N_plikow_zmierzony >= stary$N_plikow_zmierzony[[1]],
          Najlepszy_PEP <= limit_PEP,
          Max_dobrych_fragmentow >= stary$Max_dobrych_fragmentow[[1]],
          Max_sredni_score_top3 >= stary$Max_sredni_score_top3[[1]],
          Max_min_score_top3 >= stary$Max_min_score_top3[[1]],
          Max_udzial_runow_top3 >= stary$Max_udzial_runow_top3[[1]],
          Max_suma_top_fragmentow >=
            stary$Max_suma_top_fragmentow[[1]] * minimalny_udzial_intensywnosci,
          Max_srednia_top_fragmentow >=
            stary$Max_srednia_top_fragmentow[[1]] * minimalny_udzial_intensywnosci,
          Max_intensywnosc >=
            stary$Max_intensywnosc[[1]] * minimalny_udzial_intensywnosci
        )
    }
    if (nrow(alternatywy) == 0) next

    nowe_pm1 <- licz_punkty_wokol(
      alternatywy$RT_sredni,
      rt_biezace,
      polowa_okna_ruchomego_min
    )
    stary_w_nowym_pm1 <- as.integer(
      abs(alternatywy$RT_sredni - stary$RT_sredni[[1]]) <=
        polowa_okna_ruchomego_min + 1e-10
    )
    nowe_pm1 <- pmax(nowe_pm1 - stary_w_nowym_pm1 + 1L, 1L)

    nowe_okna <- floor(alternatywy$RT_sredni / szerokosc_okna_stalego_min) *
      szerokosc_okna_stalego_min
    stare_okno <- floor(stary$RT_sredni[[1]] / szerokosc_okna_stalego_min) *
      szerokosc_okna_stalego_min
    nowe_1min <- vapply(nowe_okna, function(okno) {
      obecne <- if (as.character(okno) %in% names(tab_okien)) {
        as.integer(tab_okien[[as.character(okno)]])
      } else {
        0L
      }
      obecne - as.integer(okno == stare_okno) + 1L
    }, integer(1))

    nowe_emp <- vapply(seq_len(nrow(alternatywy)), function(j) {
      n <- sum(
        starty_biezace <= alternatywy$RT_stop_sredni[[j]] + 1e-10 &
          stopy_biezace >= alternatywy$RT_start_sredni[[j]] - 1e-10
      )
      stary_naklada <- as.integer(
        stary$RT_start_sredni[[1]] <= alternatywy$RT_stop_sredni[[j]] + 1e-10 &&
          stary$RT_stop_sredni[[1]] >= alternatywy$RT_start_sredni[[j]] - 1e-10
      )
      as.integer(n - stary_naklada + 1L)
    }, integer(1))

    alternatywy <- alternatywy %>%
      mutate(
        Nowe_okno_1min = nowe_1min,
        Nowe_pm1 = nowe_pm1,
        Nowe_empiryczne = nowe_emp,
        Poprawa_okno_1min = stary$Liczba_w_oknie_1min[[1]] - Nowe_okno_1min,
        Poprawa_pm1 = stary$Liczba_wspolwystepujacych_pm1[[1]] - Nowe_pm1,
        Poprawa_empiryczna =
          stary$Liczba_nakladajacych_empirycznych[[1]] - Nowe_empiryczne
      ) %>%
      filter(
        Poprawa_okno_1min >= minimalna_poprawa_okno_1min,
        Poprawa_pm1 >= minimalna_poprawa_pm1,
        Poprawa_empiryczna >= minimalna_poprawa_empiryczna
      )
    if (nrow(alternatywy) == 0) next

    tabela_i <- alternatywy %>%
      transmute(
        Poziom_jakosci = poziom,
        Bialko,
        Nazwa_bialka,
        Stary_peptyd = stary$Peptyd[[1]],
        Nowy_peptyd = Peptyd,
        Stary_RT = stary$RT_sredni[[1]],
        Nowy_RT = RT_sredni,
        Stary_PEP = stary$Najlepszy_PEP[[1]],
        Nowy_PEP = Najlepszy_PEP,
        Mnoznik_PEP = Nowy_PEP / pmax(Stary_PEP, .Machine$double.xmin),
        Bezwzgledny_wzrost_PEP = Nowy_PEP - Stary_PEP,
        Stare_N_plikow = stary$N_plikow_zmierzony[[1]],
        Nowe_N_plikow = N_plikow_zmierzony,
        Stare_fragmenty = stary$Max_dobrych_fragmentow[[1]],
        Nowe_fragmenty = Max_dobrych_fragmentow,
        Stary_sredni_score_top3 = stary$Max_sredni_score_top3[[1]],
        Nowy_sredni_score_top3 = Max_sredni_score_top3,
        Roznica_sredniego_score_top3 = Nowy_sredni_score_top3 - Stary_sredni_score_top3,
        Stary_min_score_top3 = stary$Max_min_score_top3[[1]],
        Nowy_min_score_top3 = Max_min_score_top3,
        Roznica_min_score_top3 = Nowy_min_score_top3 - Stary_min_score_top3,
        Stary_udzial_runow_top3 = stary$Max_udzial_runow_top3[[1]],
        Nowy_udzial_runow_top3 = Max_udzial_runow_top3,
        Roznica_udzialu_runow_top3 = Nowy_udzial_runow_top3 - Stary_udzial_runow_top3,
        Udzial_sumy_fragmentow =
          Max_suma_top_fragmentow / stary$Max_suma_top_fragmentow[[1]],
        Udzial_sredniej_fragmentow =
          Max_srednia_top_fragmentow / stary$Max_srednia_top_fragmentow[[1]],
        Udzial_max_intensywnosci =
          Max_intensywnosc / stary$Max_intensywnosc[[1]],
        Stare_okno_1min = stary$Liczba_w_oknie_1min[[1]],
        Nowe_okno_1min,
        Poprawa_okno_1min,
        Stare_pm1 = stary$Liczba_wspolwystepujacych_pm1[[1]],
        Nowe_pm1,
        Poprawa_pm1,
        Stare_empiryczne = stary$Liczba_nakladajacych_empirycznych[[1]],
        Nowe_empiryczne,
        Poprawa_empiryczna,
        Nowy_pelny_ranking = Pelny_ranking_w_bialku
      )
    wyniki[[length(wyniki) + 1L]] <- tabela_i
  }

  if (length(wyniki) == 0) return(pusta_tabela_ruchow())
  bind_rows(wyniki) %>%
    arrange(
      desc(Poprawa_empiryczna),
      desc(Poprawa_pm1),
      desc(Poprawa_okno_1min),
      Mnoznik_PEP,
      desc(Roznica_sredniego_score_top3),
      desc(Roznica_min_score_top3),
      desc(Roznica_udzialu_runow_top3),
      desc(pmin(
        Udzial_sumy_fragmentow,
        Udzial_sredniej_fragmentow,
        Udzial_max_intensywnosci
      )),
      Nowy_pelny_ranking,
      Bialko,
      Stary_peptyd,
      Nowy_peptyd
    )
}

# Wykonuje jedną testową zamianę i zwraca zaktualizowaną listę.
zastap_peptyd <- function(aktualna, ruch, kandydaci) {
  nowy <- kandydaci %>%
    filter(
      Bialko == ruch$Bialko[[1]],
      Peptyd == ruch$Nowy_peptyd[[1]]
    ) %>%
    slice(1) %>%
    transmute(
      Bialko,
      Nazwa_bialka,
      Peptyd,
      Wybrany_ladunek,
      N_plikow_zmierzony,
      RT_sredni,
      RT_start_sredni,
      RT_stop_sredni,
      Najlepszy_PEP,
      Max_intensywnosc,
      Wybrane_mz,
      Max_dobrych_fragmentow,
      Max_srednia_top_fragmentow,
      Max_suma_top_fragmentow,
      Najmocniejszy_fragment,
      Max_sredni_score_top3,
      Max_min_score_top3,
      Max_udzial_runow_top3,
      Pozycja_w_bialku = as.integer(NA),
      Zakres_analizy = unique(aktualna$Zakres_analizy),
      Typ_listy = unique(aktualna$Typ_listy),
      Klucz
    )

  aktualna %>%
    filter(
      !(Bialko == ruch$Bialko[[1]] & Peptyd == ruch$Stary_peptyd[[1]])
    ) %>%
    bind_rows(nowy) %>%
    arrange(Bialko, RT_sredni, Peptyd)
}

# Iteracyjnie wybiera tylko takie zamiany, które spełniają kryteria jakości i nie pogarszają globalnego zatłoczenia.
optymalizuj_liste <- function(lista_start, kandydaci) {
  aktualna <- lista_start
  historia <- list()
  uzyte_bialka <- character(0)
  wycofane_klucze <- character(0)
  max_iter <- n_distinct(lista_start$Bialko)

  for (iter in seq_len(max_iter)) {
    ruchy <- enumeruj_ruchy(
      aktualna,
      kandydaci,
      poziom = "TOLERANCJA",
      uzyte_bialka = uzyte_bialka,
      wycofane_klucze = wycofane_klucze
    )
    if (nrow(ruchy) == 0) break

    przed <- podsumuj_zatloczenie(aktualna, "Przed ruchem")
    ocenione <- list()

    for (j in seq_len(nrow(ruchy))) {
      ruch <- ruchy[j, ]
      propozycja <- zastap_peptyd(aktualna, ruch, kandydaci)
      po <- podsumuj_zatloczenie(propozycja, "Po ruchu")

      nie_wzrastaja_maksima <-
        po$Maks_okno_1min <= przed$Maks_okno_1min &&
        po$Maks_pm1 <= przed$Maks_pm1 &&
        po$Maks_empiryczne <= przed$Maks_empiryczne

      nie_wzrastaja_sumy <-
        po$Suma_kwadratow_okien_1min <= przed$Suma_kwadratow_okien_1min &&
        po$Suma_pm1 <= przed$Suma_pm1 &&
        po$Suma_empiryczne <= przed$Suma_empiryczne

      co_najmniej_jedna_suma_spada <-
        po$Suma_kwadratow_okien_1min < przed$Suma_kwadratow_okien_1min ||
        po$Suma_pm1 < przed$Suma_pm1 ||
        po$Suma_empiryczne < przed$Suma_empiryczne

      if (!(nie_wzrastaja_maksima && nie_wzrastaja_sumy && co_najmniej_jedna_suma_spada)) {
        next
      }

      ocenione[[length(ocenione) + 1L]] <- ruch %>%
        mutate(
          Iteracja = iter,
          Maks_okno_1min_przed = przed$Maks_okno_1min,
          Maks_okno_1min_po = po$Maks_okno_1min,
          Maks_pm1_przed = przed$Maks_pm1,
          Maks_pm1_po = po$Maks_pm1,
          Maks_empiryczne_przed = przed$Maks_empiryczne,
          Maks_empiryczne_po = po$Maks_empiryczne,
          Suma_kwadratow_1min_przed = przed$Suma_kwadratow_okien_1min,
          Suma_kwadratow_1min_po = po$Suma_kwadratow_okien_1min,
          Poprawa_sumy_kwadratow_1min =
            przed$Suma_kwadratow_okien_1min - po$Suma_kwadratow_okien_1min,
          Suma_pm1_przed = przed$Suma_pm1,
          Suma_pm1_po = po$Suma_pm1,
          Poprawa_sumy_pm1 = przed$Suma_pm1 - po$Suma_pm1,
          Suma_empiryczne_przed = przed$Suma_empiryczne,
          Suma_empiryczne_po = po$Suma_empiryczne,
          Poprawa_sumy_empiryczne =
            przed$Suma_empiryczne - po$Suma_empiryczne,
          Minimalny_udzial_intensywnosci = pmin(
            Udzial_sumy_fragmentow,
            Udzial_sredniej_fragmentow,
            Udzial_max_intensywnosci
          )
        )
    }

    ocenione_df <- bind_rows(ocenione)
    if (nrow(ocenione_df) == 0) break

    najlepszy <- ocenione_df %>%
      arrange(
        desc(Poprawa_sumy_empiryczne),
        desc(Poprawa_sumy_pm1),
        desc(Poprawa_sumy_kwadratow_1min),
        desc(Poprawa_empiryczna),
        desc(Poprawa_pm1),
        desc(Poprawa_okno_1min),
        Mnoznik_PEP,
        desc(Roznica_sredniego_score_top3),
        desc(Roznica_min_score_top3),
        desc(Roznica_udzialu_runow_top3),
        desc(Minimalny_udzial_intensywnosci),
        Nowy_pelny_ranking,
        Bialko,
        Stary_peptyd,
        Nowy_peptyd
      ) %>%
      slice(1)

    aktualna <- zastap_peptyd(aktualna, najlepszy, kandydaci)
    historia[[length(historia) + 1L]] <- najlepszy
    uzyte_bialka <- c(uzyte_bialka, najlepszy$Bialko[[1]])
    wycofane_klucze <- c(
      wycofane_klucze,
      paste(najlepszy$Bialko[[1]], najlepszy$Stary_peptyd[[1]], sep = "||")
    )
  }

  historia_df <- if (length(historia) == 0) {
    pusta_tabela_ruchow() %>% mutate(Iteracja = integer())
  } else {
    bind_rows(historia)
  }

  nowe_klucze <- if (nrow(historia_df) > 0) {
    paste(historia_df$Bialko, historia_df$Nowy_peptyd, sep = "||")
  } else {
    character(0)
  }

  aktualna <- aktualna %>%
    left_join(
      kandydaci %>% select(Bialko, Peptyd, Pelny_ranking_w_bialku),
      by = c("Bialko", "Peptyd")
    ) %>%
    group_by(Bialko) %>%
    arrange(Pelny_ranking_w_bialku, Najlepszy_PEP, Peptyd, .by_group = TRUE) %>%
    mutate(
      Pozycja_w_bialku_po_propozycji = row_number(),
      Liczba_peptydow_w_bialku = n()
    ) %>%
    ungroup() %>%
    mutate(
      Status_w_propozycji = if_else(
        Klucz %in% nowe_klucze,
        "nowy kandydat — kontrolowana propozycja w tolerancji jakościowej",
        "pozostawiony z listy bazowej"
      )
    )

  list(lista = aktualna, historia = historia_df)
}

# Zapisuje wynik pojedynczej kontroli poprawności.
kontrola <- function(nazwa, wynik, obserwowana = "", oczekiwana = "") {
  tibble(
    Kontrola = nazwa,
    Wynik = isTRUE(wynik),
    Wartosc_obserwowana = as.character(obserwowana),
    Wartosc_oczekiwana = as.character(oczekiwana)
  )
}


# ------------------------------------------------------------------------------
# 4. WCZYTANIE DANYCH
# ------------------------------------------------------------------------------

# Listy bazowe pochodzą z etapu 03, natomiast pełna pula alternatywnych kandydatów z etapu 02.

sciezki_03 <- list(
  FASTA_Top2 = file.path(
    folder_etap03, "02_CSV_do_dalszego_przygotowania",
    "CSV_01A_CALA_FASTA_lista_peptydowo_prekursorowa_top2_PO_RAGGED.csv"
  ),
  FASTA_Top3 = file.path(
    folder_etap03, "02_CSV_do_dalszego_przygotowania",
    "CSV_01B_CALA_FASTA_lista_peptydowo_prekursorowa_top3_PO_RAGGED.csv"
  ),
  PANEL_Top2 = file.path(
    folder_etap03, "02_CSV_do_dalszego_przygotowania",
    "CSV_02A_PANEL_741_lista_peptydowo_prekursorowa_top2_PO_RAGGED.csv"
  ),
  PANEL_Top3 = file.path(
    folder_etap03, "02_CSV_do_dalszego_przygotowania",
    "CSV_02B_PANEL_741_lista_peptydowo_prekursorowa_top3_PO_RAGGED.csv"
  )
)

plik_kandydaci <- file.path(
  folder_etap02,
  "01_Tabele_Excel",
  "Tabela_05_PANEL_kandydaci_po_PEP_0_01_i_min_3_fragmentach.xlsx"
)

for (p in c(unlist(sciezki_03), plik_kandydaci)) {
  if (!file.exists(p)) stop("Brak wymaganego pliku: ", p)
}

listy <- list(
  FASTA_Top2 = przygotuj_liste(czytaj_csv2(sciezki_03$FASTA_Top2), "Cała FASTA", "Top 2"),
  FASTA_Top3 = przygotuj_liste(czytaj_csv2(sciezki_03$FASTA_Top3), "Cała FASTA", "Top 3"),
  PANEL_Top2 = przygotuj_liste(czytaj_csv2(sciezki_03$PANEL_Top2), "Panel 741", "Top 2"),
  PANEL_Top3 = przygotuj_liste(czytaj_csv2(sciezki_03$PANEL_Top3), "Panel 741", "Top 3")
)

kandydaci_PANEL <- przygotuj_kandydatow(readxl::read_excel(plik_kandydaci))

podsumowanie_top_etap02 <- readr::read_csv2(
  file.path(folder_etap02, "04_CSV", "ETAP_02_porownanie_Top2_Top3.csv"),
  show_col_types = FALSE,
  locale = readr::locale(decimal_mark = ",", grouping_mark = " ", encoding = "UTF-8")
) %>%
  mutate(
    Limit_ID = as.integer(na_liczbe(Limit_ID)),
    Liczba_bialek = as.integer(na_liczbe(Liczba_bialek)),
    Liczba_peptydow = as.integer(na_liczbe(Liczba_peptydow))
  )

redukcja_panel_etap02 <- readr::read_csv2(
  file.path(folder_etap02, "04_CSV", "ETAP_02_redukcja_PANEL_741.csv"),
  show_col_types = FALSE,
  locale = readr::locale(decimal_mark = ",", grouping_mark = " ", encoding = "UTF-8")
) %>%
  mutate(
    Etap_ID = as.integer(na_liczbe(Etap_ID)),
    Liczba_bialek = as.integer(na_liczbe(Liczba_bialek)),
    Liczba_peptydow = as.integer(na_liczbe(Liczba_peptydow))
  )

oczekiwane_top <- function(zbior, limit_id, kolumna) {
  x <- podsumowanie_top_etap02 %>%
    filter(Zbior == zbior, Limit_ID == limit_id) %>%
    pull(all_of(kolumna))
  if (length(x) != 1L) stop("Nie udało się odczytać wartości kontrolnej Top2/Top3 z ETAPU 02.")
  x[[1]]
}

oczekiwane_panel_fragmenty <- redukcja_panel_etap02 %>%
  filter(Etap_ID == 5L) %>%
  summarise(
    Peptydy = first(Liczba_peptydow),
    Bialka = first(Liczba_bialek)
  )

for (nm in names(listy)) {
  cat(
    nm, ": ", nrow(listy[[nm]]), " peptydów; ",
    n_distinct(listy[[nm]]$Bialko), " białek\n", sep = ""
  )
}
cat(
  "Kandydaci panelu: ", nrow(kandydaci_PANEL), " peptydów; ",
  n_distinct(kandydaci_PANEL$Bialko), " białek\n\n", sep = ""
)


# ------------------------------------------------------------------------------
# 5. ANALIZA ZATŁOCZENIA LIST BAZOWYCH
# ------------------------------------------------------------------------------

# Trzy miary zatłoczenia są liczone niezależnie, ponieważ opisują różne aspekty lokalnego nagromadzenia celów.

okna_1min <- bind_rows(lapply(listy, analiza_okien_1min))
pm1 <- bind_rows(lapply(listy, analiza_pm1))
empiryczne <- bind_rows(lapply(listy, analiza_empiryczna))
podsumowanie_bazowe <- bind_rows(lapply(listy, podsumuj_zatloczenie))

top10_okien <- okna_1min %>%
  group_by(Zakres_analizy, Typ_listy) %>%
  arrange(desc(Liczba_peptydow), Okno_start, .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

top10_pm1 <- pm1 %>%
  group_by(Zakres_analizy, Typ_listy) %>%
  arrange(desc(Liczba_wspolwystepujacych_pm1), RT_sredni, .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

top10_empiryczne <- empiryczne %>%
  group_by(Zakres_analizy, Typ_listy) %>%
  arrange(desc(Liczba_nakladajacych_empirycznych), RT_sredni, .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

bialka_1_2_PANEL <- bind_rows(
  dodaj_lokalne_zatloczenie(listy$PANEL_Top2),
  dodaj_lokalne_zatloczenie(listy$PANEL_Top3)
) %>%
  group_by(Typ_listy, Bialko) %>%
  mutate(Liczba_peptydow_w_bialku = n()) %>%
  ungroup() %>%
  filter(Liczba_peptydow_w_bialku %in% c(1L, 2L)) %>%
  arrange(Typ_listy, Liczba_peptydow_w_bialku, desc(Liczba_wspolwystepujacych_pm1))


# ------------------------------------------------------------------------------
# 6. PROPOZYCJE STRICT ORAZ KONTROLOWANE PROPOZYCJE W TOLERANCJI
# ------------------------------------------------------------------------------

# Najpierw sprawdzany jest wariant ścisły, a następnie wariant dopuszczający niewielkie różnice jakościowe zgodnie z ustalonymi progami.

strict_Top2 <- enumeruj_ruchy(listy$PANEL_Top2, kandydaci_PANEL, "STRICT")
strict_Top3 <- enumeruj_ruchy(listy$PANEL_Top3, kandydaci_PANEL, "STRICT")

wynik_Top2 <- optymalizuj_liste(listy$PANEL_Top2, kandydaci_PANEL)
wynik_Top3 <- optymalizuj_liste(listy$PANEL_Top3, kandydaci_PANEL)

panel_Top2_po <- wynik_Top2$lista
panel_Top3_po <- wynik_Top3$lista
historia_Top2 <- wynik_Top2$historia
historia_Top3 <- wynik_Top3$historia

podsumowanie_po <- bind_rows(
  podsumuj_zatloczenie(panel_Top2_po, "Po kontrolowanych propozycjach"),
  podsumuj_zatloczenie(panel_Top3_po, "Po kontrolowanych propozycjach")
)

porownanie_przed_po <- bind_rows(
  podsumowanie_bazowe %>% filter(Zakres_analizy == "Panel 741"),
  podsumowanie_po
) %>%
  arrange(Typ_listy, Wariant)


# ------------------------------------------------------------------------------
# 7. KONTROLE
# ------------------------------------------------------------------------------

# Po każdej optymalizacji weryfikowane są liczebności, zachowanie białek oraz brak pogorszenia globalnych miar zatłoczenia.

kontrole <- bind_rows(
  kontrola(
    "FASTA Top2 — liczba peptydów",
    nrow(listy$FASTA_Top2) == oczekiwane_top("Cała FASTA", 2L, "Liczba_peptydow"),
    nrow(listy$FASTA_Top2), oczekiwane_top("Cała FASTA", 2L, "Liczba_peptydow")
  ),
  kontrola(
    "FASTA Top2 — liczba białek",
    n_distinct(listy$FASTA_Top2$Bialko) == oczekiwane_top("Cała FASTA", 2L, "Liczba_bialek"),
    n_distinct(listy$FASTA_Top2$Bialko), oczekiwane_top("Cała FASTA", 2L, "Liczba_bialek")
  ),
  kontrola(
    "FASTA Top3 — liczba peptydów",
    nrow(listy$FASTA_Top3) == oczekiwane_top("Cała FASTA", 3L, "Liczba_peptydow"),
    nrow(listy$FASTA_Top3), oczekiwane_top("Cała FASTA", 3L, "Liczba_peptydow")
  ),
  kontrola(
    "FASTA Top3 — liczba białek",
    n_distinct(listy$FASTA_Top3$Bialko) == oczekiwane_top("Cała FASTA", 3L, "Liczba_bialek"),
    n_distinct(listy$FASTA_Top3$Bialko), oczekiwane_top("Cała FASTA", 3L, "Liczba_bialek")
  ),
  kontrola(
    "Panel Top2 — liczba peptydów",
    nrow(listy$PANEL_Top2) == oczekiwane_top("Panel 741 białek", 2L, "Liczba_peptydow"),
    nrow(listy$PANEL_Top2), oczekiwane_top("Panel 741 białek", 2L, "Liczba_peptydow")
  ),
  kontrola(
    "Panel Top2 — liczba białek",
    n_distinct(listy$PANEL_Top2$Bialko) == oczekiwane_top("Panel 741 białek", 2L, "Liczba_bialek"),
    n_distinct(listy$PANEL_Top2$Bialko), oczekiwane_top("Panel 741 białek", 2L, "Liczba_bialek")
  ),
  kontrola(
    "Panel Top3 — liczba peptydów",
    nrow(listy$PANEL_Top3) == oczekiwane_top("Panel 741 białek", 3L, "Liczba_peptydow"),
    nrow(listy$PANEL_Top3), oczekiwane_top("Panel 741 białek", 3L, "Liczba_peptydow")
  ),
  kontrola(
    "Panel Top3 — liczba białek",
    n_distinct(listy$PANEL_Top3$Bialko) == oczekiwane_top("Panel 741 białek", 3L, "Liczba_bialek"),
    n_distinct(listy$PANEL_Top3$Bialko), oczekiwane_top("Panel 741 białek", 3L, "Liczba_bialek")
  ),
  kontrola(
    "Pełna pula panelu — liczba kandydatów",
    nrow(kandydaci_PANEL) == oczekiwane_panel_fragmenty$Peptydy[[1]],
    nrow(kandydaci_PANEL), oczekiwane_panel_fragmenty$Peptydy[[1]]
  ),
  kontrola(
    "Pełna pula panelu — liczba białek",
    n_distinct(kandydaci_PANEL$Bialko) == oczekiwane_panel_fragmenty$Bialka[[1]],
    n_distinct(kandydaci_PANEL$Bialko), oczekiwane_panel_fragmenty$Bialka[[1]]
  ),
  kontrola(
    "Panel Top2 jest podzbiorem Panel Top3",
    all(listy$PANEL_Top2$Klucz %in% listy$PANEL_Top3$Klucz),
    sum(!listy$PANEL_Top2$Klucz %in% listy$PANEL_Top3$Klucz), 0
  ),
  kontrola(
    "FASTA Top2 jest podzbiorem FASTA Top3",
    all(listy$FASTA_Top2$Klucz %in% listy$FASTA_Top3$Klucz),
    sum(!listy$FASTA_Top2$Klucz %in% listy$FASTA_Top3$Klucz), 0
  ),
  kontrola(
    "Brak rekordu Biognosys iRT we wszystkich listach",
    !any(vapply(listy, function(x) {
      any(str_detect(paste(x$Bialko, x$Nazwa_bialka), regex("Biognosys|iRT-Kit_WR_fusion", ignore_case = TRUE)))
    }, logical(1))),
    "", "brak"
  ),
  kontrola(
    "Brak peptydów ragged we wszystkich listach",
    !any(vapply(listy, function(x) any(czy_ragged(x$Peptyd)), logical(1))),
    "", "brak"
  ),
  kontrola(
    "Brak metioniny we wszystkich listach",
    !any(vapply(listy, function(x) any(str_detect(x$Peptyd, "M")), logical(1))),
    "", "brak"
  ),
  kontrola(
    "Wszystkie listy mają PEP <= 0,01",
    all(vapply(listy, function(x) all(x$Najlepszy_PEP <= prog_PEP), logical(1))),
    "", "TRUE"
  ),
  kontrola(
    "Wszystkie listy mają >=3 kwalifikujące się fragmenty",
    all(vapply(listy, function(x) all(x$Max_dobrych_fragmentow >= min_dobrych_fragmentow), logical(1))),
    "", "TRUE"
  ),
  kontrola(
    "Top2 po propozycjach zachowuje liczbę rekordów",
    nrow(panel_Top2_po) == nrow(listy$PANEL_Top2),
    nrow(panel_Top2_po), nrow(listy$PANEL_Top2)
  ),
  kontrola(
    "Top3 po propozycjach zachowuje liczbę rekordów",
    nrow(panel_Top3_po) == nrow(listy$PANEL_Top3),
    nrow(panel_Top3_po), nrow(listy$PANEL_Top3)
  ),
  kontrola(
    "Top2 po propozycjach zachowuje zestaw białek",
    setequal(panel_Top2_po$Bialko, listy$PANEL_Top2$Bialko),
    n_distinct(panel_Top2_po$Bialko), n_distinct(listy$PANEL_Top2$Bialko)
  ),
  kontrola(
    "Top3 po propozycjach zachowuje zestaw białek",
    setequal(panel_Top3_po$Bialko, listy$PANEL_Top3$Bialko),
    n_distinct(panel_Top3_po$Bialko), n_distinct(listy$PANEL_Top3$Bialko)
  ),
  kontrola(
    "Top2 — maksymalnie jedna zamiana na białko",
    nrow(historia_Top2) == 0 || max(table(historia_Top2$Bialko)) <= maksymalna_liczba_zamian_na_bialko,
    ifelse(nrow(historia_Top2) == 0, 0, max(table(historia_Top2$Bialko))),
    maksymalna_liczba_zamian_na_bialko
  ),
  kontrola(
    "Top3 — maksymalnie jedna zamiana na białko",
    nrow(historia_Top3) == 0 || max(table(historia_Top3$Bialko)) <= maksymalna_liczba_zamian_na_bialko,
    ifelse(nrow(historia_Top3) == 0, 0, max(table(historia_Top3$Bialko))),
    maksymalna_liczba_zamian_na_bialko
  ),
  kontrola(
    "Top2 — nie zamieniono jedynego peptydu białka",
    !chron_bialka_z_jednym_peptydem || nrow(historia_Top2) == 0 ||
      all(historia_Top2$Bialko %in% (listy$PANEL_Top2 %>% count(Bialko) %>% filter(n >= 2) %>% pull(Bialko))),
    "", "TRUE"
  ),
  kontrola(
    "Top3 — nie zamieniono jedynego peptydu białka",
    !chron_bialka_z_jednym_peptydem || nrow(historia_Top3) == 0 ||
      all(historia_Top3$Bialko %in% (listy$PANEL_Top3 %>% count(Bialko) %>% filter(n >= 2) %>% pull(Bialko))),
    "", "TRUE"
  ),
  kontrola(
    "Top2 — wszystkie propozycje spełniają tolerancję jakościową",
    nrow(historia_Top2) == 0 || all(
      historia_Top2$Nowe_N_plikow >= historia_Top2$Stare_N_plikow &
        historia_Top2$Nowe_fragmenty >= historia_Top2$Stare_fragmenty &
        historia_Top2$Roznica_sredniego_score_top3 >= -1e-12 &
        historia_Top2$Roznica_min_score_top3 >= -1e-12 &
        historia_Top2$Roznica_udzialu_runow_top3 >= -1e-12 &
        historia_Top2$Mnoznik_PEP <= maksymalny_mnoznik_PEP + 1e-12 &
        historia_Top2$Bezwzgledny_wzrost_PEP <= maksymalny_bezwzgledny_wzrost_PEP + 1e-12 &
        historia_Top2$Minimalny_udzial_intensywnosci >= minimalny_udzial_intensywnosci - 1e-12
    ),
    "", "TRUE"
  ),
  kontrola(
    "Top3 — wszystkie propozycje spełniają tolerancję jakościową",
    nrow(historia_Top3) == 0 || all(
      historia_Top3$Nowe_N_plikow >= historia_Top3$Stare_N_plikow &
        historia_Top3$Nowe_fragmenty >= historia_Top3$Stare_fragmenty &
        historia_Top3$Roznica_sredniego_score_top3 >= -1e-12 &
        historia_Top3$Roznica_min_score_top3 >= -1e-12 &
        historia_Top3$Roznica_udzialu_runow_top3 >= -1e-12 &
        historia_Top3$Mnoznik_PEP <= maksymalny_mnoznik_PEP + 1e-12 &
        historia_Top3$Bezwzgledny_wzrost_PEP <= maksymalny_bezwzgledny_wzrost_PEP + 1e-12 &
        historia_Top3$Minimalny_udzial_intensywnosci >= minimalny_udzial_intensywnosci - 1e-12
    ),
    "", "TRUE"
  ),
  kontrola(
    "Top2 — każda propozycja lokalnie poprawia trzy miary zatłoczenia",
    nrow(historia_Top2) == 0 || all(
      historia_Top2$Poprawa_okno_1min >= minimalna_poprawa_okno_1min &
        historia_Top2$Poprawa_pm1 >= minimalna_poprawa_pm1 &
        historia_Top2$Poprawa_empiryczna >= minimalna_poprawa_empiryczna
    ),
    "", "TRUE"
  ),
  kontrola(
    "Top3 — każda propozycja lokalnie poprawia trzy miary zatłoczenia",
    nrow(historia_Top3) == 0 || all(
      historia_Top3$Poprawa_okno_1min >= minimalna_poprawa_okno_1min &
        historia_Top3$Poprawa_pm1 >= minimalna_poprawa_pm1 &
        historia_Top3$Poprawa_empiryczna >= minimalna_poprawa_empiryczna
    ),
    "", "TRUE"
  )
)

# Globalne porównanie przed/po.
przed_Top2 <- podsumowanie_bazowe %>% filter(Zakres_analizy == "Panel 741", Typ_listy == "Top 2")
przed_Top3 <- podsumowanie_bazowe %>% filter(Zakres_analizy == "Panel 741", Typ_listy == "Top 3")
po_Top2 <- podsumowanie_po %>% filter(Typ_listy == "Top 2")
po_Top3 <- podsumowanie_po %>% filter(Typ_listy == "Top 3")

kontrole <- bind_rows(
  kontrole,
  kontrola(
    "Top2 — globalne maksima zatłoczenia nie wzrosły",
    po_Top2$Maks_okno_1min <= przed_Top2$Maks_okno_1min &&
      po_Top2$Maks_pm1 <= przed_Top2$Maks_pm1 &&
      po_Top2$Maks_empiryczne <= przed_Top2$Maks_empiryczne,
    paste(po_Top2$Maks_okno_1min, po_Top2$Maks_pm1, po_Top2$Maks_empiryczne, sep = "/"),
    paste(przed_Top2$Maks_okno_1min, przed_Top2$Maks_pm1, przed_Top2$Maks_empiryczne, sep = "/")
  ),
  kontrola(
    "Top3 — globalne maksima zatłoczenia nie wzrosły",
    po_Top3$Maks_okno_1min <= przed_Top3$Maks_okno_1min &&
      po_Top3$Maks_pm1 <= przed_Top3$Maks_pm1 &&
      po_Top3$Maks_empiryczne <= przed_Top3$Maks_empiryczne,
    paste(po_Top3$Maks_okno_1min, po_Top3$Maks_pm1, po_Top3$Maks_empiryczne, sep = "/"),
    paste(przed_Top3$Maks_okno_1min, przed_Top3$Maks_pm1, przed_Top3$Maks_empiryczne, sep = "/")
  ),
  kontrola(
    "Top2 — globalne sumy zatłoczenia nie wzrosły",
    po_Top2$Suma_kwadratow_okien_1min <= przed_Top2$Suma_kwadratow_okien_1min &&
      po_Top2$Suma_pm1 <= przed_Top2$Suma_pm1 &&
      po_Top2$Suma_empiryczne <= przed_Top2$Suma_empiryczne,
    paste(po_Top2$Suma_kwadratow_okien_1min, po_Top2$Suma_pm1, po_Top2$Suma_empiryczne, sep = "/"),
    paste(przed_Top2$Suma_kwadratow_okien_1min, przed_Top2$Suma_pm1, przed_Top2$Suma_empiryczne, sep = "/")
  ),
  kontrola(
    "Top3 — globalne sumy zatłoczenia nie wzrosły",
    po_Top3$Suma_kwadratow_okien_1min <= przed_Top3$Suma_kwadratow_okien_1min &&
      po_Top3$Suma_pm1 <= przed_Top3$Suma_pm1 &&
      po_Top3$Suma_empiryczne <= przed_Top3$Suma_empiryczne,
    paste(po_Top3$Suma_kwadratow_okien_1min, po_Top3$Suma_pm1, po_Top3$Suma_empiryczne, sep = "/"),
    paste(przed_Top3$Suma_kwadratow_okien_1min, przed_Top3$Suma_pm1, przed_Top3$Suma_empiryczne, sep = "/")
  )
)

if (any(!kontrole$Wynik)) {
  print(kontrole %>% filter(!Wynik))
  stop("ETAP 04 zatrzymany: co najmniej jedna kontrola krytyczna nie została spełniona.")
}


# ------------------------------------------------------------------------------
# 8. EKSPORT TABEL
# ------------------------------------------------------------------------------

zapisz_wiele_arkuszy(
  list(
    Podsumowanie = bind_rows(podsumowanie_bazowe, podsumowanie_po),
    Top10_okna_1min = top10_okien,
    Top10_pm1 = top10_pm1,
    Top10_empiryczne = top10_empiryczne
  ),
  file.path(folder_tabele, "Tabela_01_Podsumowanie_zatloczenia_Top2_Top3.xlsx")
)

zapisz_wiele_arkuszy(
  list(
    Okna_1min = okna_1min,
    Otoczenie_pm1 = pm1,
    Empiryczne_RT = empiryczne
  ),
  file.path(folder_tabele, "Tabela_02_Pelna_analiza_okien_RETENCJI.xlsx")
)

zapisz_excel(
  bialka_1_2_PANEL,
  file.path(folder_tabele, "Tabela_03_PANEL_bialka_z_1_lub_2_peptydami.xlsx"),
  "Bialka_1_2"
)

zapisz_wiele_arkuszy(
  list(
    STRICT_Top2 = strict_Top2,
    STRICT_Top3 = strict_Top3
  ),
  file.path(folder_tabele, "Tabela_04_Strict_alternatywy_bez_pogorszenia_jakosci.xlsx")
)

zapisz_wiele_arkuszy(
  list(
    Historia_Top2 = historia_Top2,
    Historia_Top3 = historia_Top3
  ),
  file.path(folder_tabele, "Tabela_05_Historia_kontrolowanych_propozycji_zamian.xlsx")
)

zapisz_wiele_arkuszy(
  list(
    Panel_Top2_po = panel_Top2_po,
    Panel_Top3_po = panel_Top3_po
  ),
  file.path(folder_tabele, "Tabela_06_PANEL_listy_po_kontrolowanych_propozycjach.xlsx")
)

zapisz_excel(
  porownanie_przed_po,
  file.path(folder_tabele, "Tabela_07_PANEL_porownanie_zatloczenia_przed_po.xlsx"),
  "Przed_po"
)

zapisz_csv2(okna_1min, "ETAP_04_okna_1min_Top2_Top3.csv")
zapisz_csv2(pm1, "ETAP_04_otoczenie_RT_pm1_Top2_Top3.csv")
zapisz_csv2(empiryczne, "ETAP_04_nakladanie_empirycznych_okien_Top2_Top3.csv")
zapisz_csv2(bind_rows(podsumowanie_bazowe, podsumowanie_po), "ETAP_04_podsumowanie_zatloczenia.csv")
zapisz_csv2(strict_Top2, "ETAP_04_PANEL_Top2_alternatywy_STRICT.csv")
zapisz_csv2(strict_Top3, "ETAP_04_PANEL_Top3_alternatywy_STRICT.csv")
zapisz_csv2(historia_Top2, "ETAP_04_PANEL_Top2_historia_kontrolowanych_zamian.csv")
zapisz_csv2(historia_Top3, "ETAP_04_PANEL_Top3_historia_kontrolowanych_zamian.csv")
zapisz_csv2(panel_Top2_po, "ETAP_04_PANEL_Top2_lista_po_kontrolowanych_propozycjach.csv")
zapisz_csv2(panel_Top3_po, "ETAP_04_PANEL_Top3_lista_po_kontrolowanych_propozycjach.csv")

# Preferowany alias odpowiadający panelowej liście Top 3 po propozycjach zamian.
zapisz_csv2(
  panel_Top3_po,
  "ETAP_04_PANEL_lista_po_propozycjach_zamian.csv"
)

# Historyczny alias zgodnościowy dla starszych skryptów.
zapisz_csv2(
  panel_Top3_po,
  "ETAP_04_PANEL_lista_po_eksploracyjnej_propozycji.csv"
)

zapisz_excel(
  kontrole,
  file.path(folder_kontrola, "KONTROLA_ETAP_04_FINAL.xlsx"),
  "Kontrole"
)
readr::write_excel_csv2(
  kontrole,
  file.path(folder_kontrola, "KONTROLA_ETAP_04_FINAL.csv")
)


# ------------------------------------------------------------------------------
# 9. WYKRESY
# ------------------------------------------------------------------------------

theme_prm <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2),
      axis.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(15, 25, 20, 20)
    )
}

zapisz_wykres <- function(p, nazwa, szer = 11, wys = 6.5) {
  ggplot2::ggsave(
    file.path(folder_wykresy, paste0(nazwa, ".png")),
    p, width = szer, height = wys, dpi = 300, bg = "white"
  )
  ggplot2::ggsave(
    file.path(folder_wykresy, paste0(nazwa, ".pdf")),
    p, width = szer, height = wys, bg = "white"
  )
}

p1 <- ggplot(okna_1min, aes(x = Okno_start, y = Liczba_peptydow)) +
  geom_col() +
  facet_grid(Zakres_analizy ~ Typ_listy, scales = "free_y") +
  labs(
    title = "Rozmieszczenie peptydów w 1-minutowych przedziałach RT",
    x = "Początek przedziału RT [min]",
    y = "Liczba peptydów"
  ) +
  theme_prm()
zapisz_wykres(p1, "Wykres_01_Okna_1min_FASTA_PANEL_Top2_Top3", 12, 8)

p2 <- ggplot(
  porownanie_przed_po,
  aes(x = Typ_listy, y = Maks_pm1, fill = Wariant)
) +
  geom_col(position = "dodge") +
  labs(
    title = "Panel: maksymalne zatłoczenie w otoczeniu RT ±1 min",
    x = NULL,
    y = "Maksymalna liczba współwystępujących peptydów",
    fill = "Wariant"
  ) +
  theme_prm()
zapisz_wykres(p2, "Wykres_02_PANEL_maks_pm1_przed_po", 9, 6)

p3 <- ggplot(
  porownanie_przed_po,
  aes(x = Typ_listy, y = Maks_empiryczne, fill = Wariant)
) +
  geom_col(position = "dodge") +
  labs(
    title = "Panel: maksymalne nakładanie empirycznych przedziałów RT",
    x = NULL,
    y = "Maksymalna liczba nakładających się przedziałów",
    fill = "Wariant"
  ) +
  theme_prm()
zapisz_wykres(p3, "Wykres_03_PANEL_empiryczne_przed_po", 9, 6)

liczby_zamian <- tibble(
  Typ_listy = c("Top 2", "Top 3"),
  Propozycje_STRICT = c(nrow(strict_Top2), nrow(strict_Top3)),
  Zastosowane_w_tolerancji = c(nrow(historia_Top2), nrow(historia_Top3))
) %>%
  pivot_longer(-Typ_listy, names_to = "Rodzaj", values_to = "Liczba")

p4 <- ggplot(liczby_zamian, aes(x = Typ_listy, y = Liczba, fill = Rodzaj)) +
  geom_col(position = "dodge") +
  labs(
    title = "Liczba jakościowo dopuszczalnych propozycji zamiany",
    x = NULL,
    y = "Liczba propozycji",
    fill = "Rodzaj"
  ) +
  theme_prm()
zapisz_wykres(p4, "Wykres_04_Liczba_propozycji_zamian", 9, 6)


# ------------------------------------------------------------------------------
# 10. INDEKS I RAPORT
# ------------------------------------------------------------------------------

indeks <- tribble(
  ~Folder, ~Plik, ~Opis,
  "01_Tabele_Excel", "Tabela_01_Podsumowanie_zatloczenia_Top2_Top3.xlsx", "Podsumowanie i najbardziej zatłoczone zakresy dla FASTA/panelu i Top2/Top3.",
  "01_Tabele_Excel", "Tabela_02_Pelna_analiza_okien_RETENCJI.xlsx", "Pełne wyniki okien 1 min, RT ±1 min i nakładania empirycznego.",
  "01_Tabele_Excel", "Tabela_03_PANEL_bialka_z_1_lub_2_peptydami.xlsx", "Kontrola białek panelowych reprezentowanych jednym lub dwoma peptydami.",
  "01_Tabele_Excel", "Tabela_04_Strict_alternatywy_bez_pogorszenia_jakosci.xlsx", "Alternatywy, które nie pogarszają żadnej kontrolowanej miary jakości.",
  "01_Tabele_Excel", "Tabela_05_Historia_kontrolowanych_propozycji_zamian.xlsx", "Kolejne ruchy dopuszczone w ściśle zdefiniowanej tolerancji jakościowej.",
  "01_Tabele_Excel", "Tabela_06_PANEL_listy_po_kontrolowanych_propozycjach.xlsx", "Warianty Top2 i Top3 po kontrolowanych propozycjach.",
  "01_Tabele_Excel", "Tabela_07_PANEL_porownanie_zatloczenia_przed_po.xlsx", "Porównanie globalnych miar zatłoczenia przed i po.",
  "03_Kontrola", "KONTROLA_ETAP_04_FINAL.xlsx/csv", "Kontrole liczebności, jakości, integralności i braku pogorszenia zatłoczenia.",
  "04_CSV", "ETAP_04_PANEL_Top2_lista_po_kontrolowanych_propozycjach.csv", "Wariant panelowej listy Top2.",
  "04_CSV", "ETAP_04_PANEL_Top3_lista_po_kontrolowanych_propozycjach.csv", "Wariant panelowej listy Top3."
)

zapisz_excel(
  indeks,
  file.path(folder_wyniki, "00_INDEKS_PLIKOW_ETAP_04_FINAL.xlsx"),
  "Indeks"
)

raport_path <- file.path(folder_kontrola, "RAPORT_ETAP_04_FINAL.txt")
sink(raport_path)
cat("ETAP 04 — OKNA RETENCJI I KONTROLOWANE PROPOZYCJE ZAMIANY\n")
cat("Data wykonania:", format(Sys.time()), "\n\n")
cat("Źródło ETAP 03:", folder_etap03, "\n")
cat("Źródło ETAP 02:", folder_etap02, "\n\n")
cat("PARAMETRY TOLERANCJI JAKOŚCIOWEJ\n")
cat("Maksymalny mnożnik PEP:", maksymalny_mnoznik_PEP, "\n")
cat("Maksymalny bezwzględny wzrost PEP:", maksymalny_bezwzgledny_wzrost_PEP, "\n")
cat("Minimalny udział każdej miary intensywności:", minimalny_udzial_intensywnosci, "\n")
cat("Liczba plików: nie mniejsza\n")
cat("Liczba kwalifikujących się fragmentów: nie mniejsza\n")
cat("Średni Fr.Score top 3: nie mniejszy\n")
cat("Minimalny Fr.Score top 3: nie mniejszy\n")
cat("Udział runów top 3: nie mniejszy\n")
cat("Chronienie białek z jednym peptydem:", chron_bialka_z_jednym_peptydem, "\n\n")
cat("PROPOZYCJE STRICT\n")
cat("Top 2:", nrow(strict_Top2), "\n")
cat("Top 3:", nrow(strict_Top3), "\n\n")
cat("KONTROLOWANE PROPOZYCJE W TOLERANCJI\n")
cat("Top 2:", nrow(historia_Top2), "\n")
cat("Top 3:", nrow(historia_Top3), "\n\n")
cat("PODSUMOWANIE ZATŁOCZENIA\n")
print(bind_rows(podsumowanie_bazowe, podsumowanie_po))
cat("\nWszystkie kontrole krytyczne: TRUE\n\n")
cat(
  "Interpretacja: warianty po zamianach są propozycjami do dalszej oceny. ",
  "Każda zastosowana zamiana zachowuje liczbę plików i kwalifikujących się fragmentów, ",
  "nie obniża metryk Fr.Score/powtarzalności top 3, ",
  "spełnia jawnie zdefiniowaną tolerancję PEP i intensywności oraz ",
  "nie zwiększa żadnej globalnej miary zatłoczenia. ",
  "Nie należy traktować ich jako automatycznie zatwierdzonej listy PRM.\n",
  sep = ""
)
sink()

cat("\nGOTOWE — ETAP 04\n")
cat("Folder wynikowy:\n", normalizePath(folder_wyniki, winslash = "/"), "\n")
cat("Propozycje STRICT Top2/Top3:", nrow(strict_Top2), "/", nrow(strict_Top3), "\n")
cat("Kontrolowane propozycje Top2/Top3:", nrow(historia_Top2), "/", nrow(historia_Top3), "\n")

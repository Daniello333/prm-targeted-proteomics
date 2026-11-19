# ==============================================================================
# 05. KALIBRACJA iRT I LISTY ROZSZERZONE
#
# Praca magisterska:
# „Planowanie celowanych eksperymentów proteomicznych z wykorzystaniem środowiska R”
# Autor: Daniel Drabik
#
# Skrypt wykorzystuje dane 30- i 60-minutowe do oceny transferu czasu retencji
# przez skalę iRT. Bazowe listy Top 2 i Top 3 pochodzą z etapu 03. Wolne miejsca
# mogą być uzupełniane kandydatami spełniającymi kryteria jakości w danych
# 60-minutowych.
#
# Jeżeli peptyd był obserwowany w którymkolwiek pomiarze 30-minutowym, zachowany
# zostaje jego obserwowany RT30. Przewidywany RT30 jest wyznaczany tylko dla
# peptydów niezaobserwowanych w żadnym pliku 30-minutowym. Otrzymane listy są
# przeznaczone do dalszej kontroli i planowania metody PRM.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. PAKIETY
# ------------------------------------------------------------------------------

pakiety <- c("tidyverse", "openxlsx", "ggplot2")

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

prog_PEP <- 0.01
min_dobrych_fragmentow <- 3L
dozwolone_ladunki <- c(2L, 3L)

# Liczebności wyników nie są wpisane na sztywno.
# Kontrole wejścia są odczytywane z finalnych podsumowań ETAPU 02,
# a liczebności nowych/rozszerzonych list są zapisywane jako dynamiczny snapshot.

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
  paste0(
    "WYNIKI_ETAP_05_KALKULATOR_iRT_I_LISTY_ROZSZERZONE_EXPORT_QUANT_HUMAN_ONLY_",
    znacznik
  )
)

folder_tabele <- file.path(folder_wyniki, "01_Tabele_Excel")
folder_csv <- file.path(folder_wyniki, "02_CSV")
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
# 2. WYSZUKIWANIE FOLDERÓW LUB ZIP-ÓW WEJŚCIOWYCH
# ------------------------------------------------------------------------------

# Sprawdza komplet plików wymaganych z wcześniejszych etapów.
ma_pliki <- function(folder, wymagane) {
  dir.exists(folder) && all(file.exists(file.path(folder, wymagane)))
}

# Wybiera najnowszy zgodny katalog lub archiwum.
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
  katalog_temp <- tempfile(pattern = paste0("ETAP05_", etykieta, "_"))
  dir.create(katalog_temp, recursive = TRUE, showWarnings = FALSE)
  utils::unzip(zip_path, exdir = katalog_temp)

  kandydaci <- c(
    katalog_temp,
    list.dirs(katalog_temp, recursive = TRUE, full.names = TRUE)
  )

  kandydaci <- kandydaci[
    vapply(kandydaci, ma_pliki, logical(1), wymagane = wymagane)
  ]

  wybierz_najnowszy(kandydaci)
}

# Odnajduje komplet danych z etapu 02 lub 03, także wewnątrz archiwum ZIP.
znajdz_lub_wypakuj <- function(
  wymagane,
  preferowany_wzorzec_folderu,
  preferowany_wzorzec_zip,
  etykieta
) {
  foldery <- lista_folderow()
  kandydaci <- foldery[
    vapply(foldery, ma_pliki, logical(1), wymagane = wymagane)
  ]

  preferowane <- kandydaci[stringr::str_detect(
    basename(kandydaci),
    stringr::regex(preferowany_wzorzec_folderu, ignore_case = TRUE)
  )]

  if (length(preferowane) > 0) return(wybierz_najnowszy(preferowane))
  if (length(kandydaci) > 0) return(wybierz_najnowszy(kandydaci))

  zipy <- lista_zipow()
  zipy <- zipy[stringr::str_detect(
    basename(zipy),
    stringr::regex(preferowany_wzorzec_zip, ignore_case = TRUE)
  )]

  zipy <- zipy[
    vapply(zipy, zip_zawiera, logical(1), wymagane = wymagane)
  ]

  if (length(zipy) == 0) return(NULL)

  wypakuj_do_temp(
    wybierz_najnowszy(zipy),
    wymagane,
    etykieta
  )
}

wymagane_etap03 <- c(
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
    "04_CSV",
    "ETAP_02_status_741_bialek_panelu.csv"
  ),
  file.path(
    "04_CSV",
    "ETAP_02_porownanie_Top2_Top3.csv"
  ),
  file.path(
    "04_CSV",
    "ETAP_02_porownanie_0_vs_1_DPMC_30_60_szczegolowe.csv"
  ),
  file.path(
    "00_DANE_POSREDNIE_DIANN",
    "E2_Trypsyna_1MC_RT_30min.csv"
  ),
  file.path(
    "00_DANE_POSREDNIE_DIANN",
    "E2_Trypsyna_1MC_RT_60min.csv"
  ),
  file.path(
    "00_DANE_POSREDNIE_DIANN",
    "E4_kalibracja_iRT.csv"
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
cat("E2 i kalibracja iRT z finalnego ETAPU 02:\n", folder_etap02, "\n\n")
cat("Wyniki ETAPU 05:\n", folder_wyniki, "\n\n")


# ------------------------------------------------------------------------------
# 3. FUNKCJE POMOCNICZE
# ------------------------------------------------------------------------------

# Ujednolica wartości liczbowe pochodzące z plików CSV.
na_liczbe <- function(x) {
  if (is.numeric(x)) return(x)

  readr::parse_number(
    as.character(x),
    locale = readr::locale(
      decimal_mark = ",",
      grouping_mark = " "
    )
  )
}

# Konwertuje wartości logiczne zapisane jako tekst.
na_logiczna <- function(x) {
  if (is.logical(x)) return(x)
  toupper(as.character(x)) %in% c("TRUE", "T", "1", "TAK")
}

czy_ragged <- function(peptyd) {
  stringr::str_detect(peptyd, "^(KK|KR|RK|RR)") |
    stringr::str_detect(peptyd, "(KK|KR|RK|RR)$")
}

# Rozpoznaje rekord techniczny standardów iRT, aby nie traktować go jako białka panelowego.
czy_standard_iRT <- function(bialko) {
  stringr::str_detect(
    as.character(bialko),
    stringr::regex(
      "^Biognosys[|]iRT-Kit_WR_fusion$",
      ignore_case = TRUE
    )
  )
}

# Sprawdza obecność kolumn wymaganych w kolejnych obliczeniach.
sprawdz_kolumny <- function(df, wymagane, nazwa) {
  brak <- setdiff(wymagane, colnames(df))
  if (length(brak) > 0) {
    stop(
      "W tabeli ", nazwa, " brakuje kolumn: ",
      paste(brak, collapse = ", ")
    )
  }
}

zapisz_excel <- function(df, sciezka) {
  openxlsx::write.xlsx(df, file = sciezka, overwrite = TRUE)
}

zapisz_excel_wiele <- function(lista_arkuszy, sciezka) {
  wb <- openxlsx::createWorkbook()

  for (nm in names(lista_arkuszy)) {
    arkusz <- substr(nm, 1, 31)
    dane <- lista_arkuszy[[nm]]

    openxlsx::addWorksheet(wb, arkusz)
    openxlsx::writeData(wb, arkusz, dane)
    openxlsx::freezePane(wb, arkusz, firstRow = TRUE)

    if (ncol(dane) > 0) {
      openxlsx::setColWidths(
        wb,
        arkusz,
        cols = seq_len(ncol(dane)),
        widths = "auto"
      )
    }
  }

  openxlsx::saveWorkbook(wb, sciezka, overwrite = TRUE)
}

zapisz_csv2 <- function(df, nazwa) {
  readr::write_excel_csv2(
    df,
    file.path(folder_csv, nazwa)
  )
}

# Dodaje wynik testu do tabeli kontroli poprawności.
kontrola <- function(nazwa, wynik, obserwowana, oczekiwana, krytyczna = TRUE) {
  tibble(
    Kontrola = nazwa,
    Wynik = isTRUE(wynik),
    Wartosc_obserwowana = as.character(obserwowana),
    Wartosc_oczekiwana = as.character(oczekiwana),
    Krytyczna = krytyczna
  )
}

theme_prm <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = base_size + 3
      ),
      plot.subtitle = ggplot2::element_text(size = base_size),
      axis.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.caption = ggplot2::element_text(
        size = base_size - 2,
        hjust = 0
      )
    )
}

zapisz_wykres <- function(wykres, nazwa, szerokosc = 10, wysokosc = 6) {
  ggplot2::ggsave(
    file.path(folder_wykresy, paste0(nazwa, ".png")),
    plot = wykres,
    width = szerokosc,
    height = wysokosc,
    dpi = 300,
    bg = "white"
  )

  ggplot2::ggsave(
    file.path(folder_wykresy, paste0(nazwa, ".pdf")),
    plot = wykres,
    width = szerokosc,
    height = wysokosc,
    bg = "white"
  )
}

# Porządkuje kandydatów według tych samych cech jakościowych, które są wykorzystywane w końcowej selekcji.
ranking_kandydatow <- function(df) {
  df %>%
    mutate(
      .rank_pliki = replace_na(N_plikow_zmierzony, 0L),
      .rank_PEP = replace_na(Najlepszy_PEP, Inf),
      .rank_score = replace_na(Max_sredni_score_top3, -Inf),
      .rank_suma = replace_na(Max_suma_top_fragmentow, 0),
      .rank_srednia = replace_na(Max_srednia_top_fragmentow, 0),
      .rank_int = replace_na(Max_intensywnosc, 0),
      .rank_fragmenty = replace_na(Max_dobrych_fragmentow, 0L)
    ) %>%
    group_by(Bialko) %>%
    arrange(
      desc(.rank_pliki),
      .rank_PEP,
      desc(.rank_score),
      desc(.rank_suma),
      desc(.rank_srednia),
      desc(.rank_int),
      desc(.rank_fragmenty),
      Peptyd,
      .by_group = TRUE
    ) %>%
    mutate(Ranking_w_puli = row_number()) %>%
    ungroup() %>%
    select(-starts_with(".rank_"))
}


# ------------------------------------------------------------------------------
# 4. WCZYTANIE KALIBRACJI iRT
# ------------------------------------------------------------------------------

# Modele kalibracyjne są odczytywane z danych wygenerowanych w etapie 02.

plik_kalibracja <- file.path(folder_etap02, "00_DANE_POSREDNIE_DIANN", "E4_kalibracja_iRT.csv")

kalibracja <- readr::read_csv2(
  plik_kalibracja,
  show_col_types = FALSE,
  locale = readr::locale(
    decimal_mark = ",",
    grouping_mark = " ",
    encoding = "UTF-8"
  )
) %>%
  mutate(
    Gradient = as.character(Gradient),
    Nachylenie_a = na_liczbe(Nachylenie_a),
    Przeciecie_b = na_liczbe(Przeciecie_b),
    R2 = na_liczbe(R2),
    N_standardow = as.integer(na_liczbe(N_standardow))
  )

sprawdz_kolumny(
  kalibracja,
  c("Gradient", "Nachylenie_a", "Przeciecie_b", "R2", "N_standardow"),
  "E4_kalibracja_iRT.csv"
)

if (!all(c("30min", "60min") %in% kalibracja$Gradient)) {
  stop("Kalibracja musi zawierać wiersze 30min i 60min.")
}

a30 <- kalibracja$Nachylenie_a[kalibracja$Gradient == "30min"][1]
b30 <- kalibracja$Przeciecie_b[kalibracja$Gradient == "30min"][1]
a60 <- kalibracja$Nachylenie_a[kalibracja$Gradient == "60min"][1]
b60 <- kalibracja$Przeciecie_b[kalibracja$Gradient == "60min"][1]

# Przelicza obserwowany RT30 na skalę iRT na podstawie modelu kalibracyjnego.
RT30_na_iRT <- function(rt30) {
  a30 * rt30 + b30
}

# Przelicza obserwowany RT60 na skalę iRT.
RT60_na_iRT <- function(rt60) {
  a60 * rt60 + b60
}

# Przelicza wartość iRT na przewidywany czas retencji w gradiencie 30-minutowym.
iRT_na_RT30 <- function(irt) {
  (irt - b30) / a30
}

# Przenosi RT60 do warunków gradientu 30-minutowego przez wspólną skalę iRT.
RT60_na_RT30 <- function(rt60) {
  iRT_na_RT30(RT60_na_iRT(rt60))
}

kalibracja_opis <- kalibracja %>%
  mutate(
    Rownanie = paste0(
      "iRT = ",
      round(Nachylenie_a, 8),
      " * RT + ",
      round(Przeciecie_b, 8)
    ),
    Pochodzenie = paste(
      "Parametry z finalnego ETAPU 02 (00_DANE_POSREDNIE_DIANN/E4_kalibracja_iRT.csv),",
      "wyznaczone na podstawie standardów iRT w danych DIA-NN."
    )
  )

wzory <- tibble(
  Element = c(
    "Model 30 min",
    "Model 60 min",
    "RT60 -> iRT",
    "iRT -> RT30",
    "RT60 -> RT30"
  ),
  Wzor = c(
    "iRT = a30 * RT30 + b30",
    "iRT = a60 * RT60 + b60",
    "iRT = a60 * RT60 + b60",
    "RT30 = (iRT - b30) / a30",
    "RT30_pred = (a60 * RT60 + b60 - b30) / a30"
  ),
  Interpretacja = c(
    "Kalibracja gradientu 30-minutowego.",
    "Kalibracja gradientu 60-minutowego.",
    "Przeniesienie obserwowanego RT60 na wspólną skalę iRT.",
    "Przeniesienie wartości iRT do gradientu 30-minutowego.",
    "Predykcja nie jest prostym podzieleniem RT60 przez dwa."
  )
)


# ------------------------------------------------------------------------------
# 5. WCZYTANIE E2 DLA 30 I 60 MIN
# ------------------------------------------------------------------------------

# Tabele obu gradientów są przygotowywane w tej samej postaci, aby można było bezpośrednio porównywać obserwacje.

kolumny_E2 <- c(
  "Bialko",
  "Nazwa_bialka",
  "Peptyd",
  "Wybrany_ladunek",
  "N_plikow_zmierzony",
  "RT_sredni",
  "RT_start_sredni",
  "RT_stop_sredni",
  "Najlepszy_PEP",
  "Max_intensywnosc",
  "Wybrane_mz",
  "Max_dobrych_fragmentow",
  "Max_srednia_top_fragmentow",
  "Max_suma_top_fragmentow",
  "Najmocniejszy_fragment",
  "Max_sredni_score_top3",
  "Max_min_score_top3",
  "Max_udzial_runow_top3",
  "Zmierzony"
)

# Wczytuje tabelę E2 dla danego gradientu i przygotowuje kolumny potrzebne do kwalifikacji kandydatów.
wczytaj_E2 <- function(sciezka, gradient) {
  cat("Wczytuję E2", gradient, "...\n")

  df <- readr::read_csv2(
    sciezka,
    show_col_types = FALSE,
    locale = readr::locale(
      decimal_mark = ",",
      grouping_mark = " ",
      encoding = "UTF-8"
    ),
    col_select = dplyr::any_of(kolumny_E2),
    na = c("", "NA", "NaN")
  )

  sprawdz_kolumny(df, kolumny_E2, paste0("E2 ", gradient))

  df %>%
    mutate(
      Bialko = as.character(Bialko),
      Nazwa_bialka = as.character(Nazwa_bialka),
      Peptyd = as.character(Peptyd),
      Wybrany_ladunek = as.integer(na_liczbe(Wybrany_ladunek)),
      N_plikow_zmierzony = as.integer(na_liczbe(N_plikow_zmierzony)),
      RT_sredni = na_liczbe(RT_sredni),
      RT_start_sredni = na_liczbe(RT_start_sredni),
      RT_stop_sredni = na_liczbe(RT_stop_sredni),
      Najlepszy_PEP = na_liczbe(Najlepszy_PEP),
      Max_intensywnosc = na_liczbe(Max_intensywnosc),
      Wybrane_mz = na_liczbe(Wybrane_mz),
      Max_dobrych_fragmentow = as.integer(
        na_liczbe(Max_dobrych_fragmentow)
      ),
      Max_srednia_top_fragmentow = na_liczbe(
        Max_srednia_top_fragmentow
      ),
      Max_suma_top_fragmentow = na_liczbe(
        Max_suma_top_fragmentow
      ),
      Najmocniejszy_fragment = na_liczbe(Najmocniejszy_fragment),
      Max_sredni_score_top3 = na_liczbe(Max_sredni_score_top3),
      Max_min_score_top3 = na_liczbe(Max_min_score_top3),
      Max_udzial_runow_top3 = na_liczbe(Max_udzial_runow_top3),
      Zmierzony = na_logiczna(Zmierzony),
      Gradient = gradient,
      Czy_ragged = czy_ragged(Peptyd),
      Czy_zawiera_M = stringr::str_detect(Peptyd, "M"),
      Czy_standard_iRT = czy_standard_iRT(Bialko)
    ) %>%
    distinct(Bialko, Peptyd, .keep_all = TRUE)
}

e2_30_raw <- wczytaj_E2(
  file.path(folder_etap02, "00_DANE_POSREDNIE_DIANN", "E2_Trypsyna_1MC_RT_30min.csv"),
  "30min"
)

e2_60_raw <- wczytaj_E2(
  file.path(folder_etap02, "00_DANE_POSREDNIE_DIANN", "E2_Trypsyna_1MC_RT_60min.csv"),
  "60min"
)

# Pozostawia rzeczywiście zaobserwowane prekursory.
kwalifikuj_obserwowane <- function(df) {
  df %>%
    filter(
      Zmierzony,
      !is.na(Bialko),
      !is.na(Peptyd)
    )
}

# Wymaga poprawnego czasu retencji oraz stanu naładowania 2+ lub 3+.
kwalifikuj_RT_i_ladunek <- function(df) {
  df %>%
    filter(
      !Czy_standard_iRT,
      !Czy_ragged,
      !Czy_zawiera_M,
      Wybrany_ladunek %in% dozwolone_ladunki,
      is.finite(RT_sredni),
      is.finite(RT_start_sredni),
      is.finite(RT_stop_sredni)
    )
}

# Stosuje kryterium PEP i minimalną liczbę kwalifikujących się fragmentów.
filtr_jakosci <- function(df) {
  df %>%
    filter(
      is.finite(Najlepszy_PEP),
      Najlepszy_PEP <= prog_PEP,
      !is.na(Max_dobrych_fragmentow),
      Max_dobrych_fragmentow >= min_dobrych_fragmentow
    )
}

obserwowane_30_all <- kwalifikuj_obserwowane(e2_30_raw)
obserwowane_60_all <- kwalifikuj_obserwowane(e2_60_raw)

kwalifikowane_30_all <- obserwowane_30_all %>% kwalifikuj_RT_i_ladunek()
kwalifikowane_60_all <- obserwowane_60_all %>% kwalifikuj_RT_i_ladunek()

jakosc_30_all <- kwalifikowane_30_all %>% filtr_jakosci()
jakosc_60_all <- kwalifikowane_60_all %>% filtr_jakosci()


# ------------------------------------------------------------------------------
# 6. PANEL 741 I LISTY BAZOWE TOP 2 / TOP 3
# ------------------------------------------------------------------------------

status_panel <- readr::read_csv2(
  file.path(
    folder_etap02,
    "04_CSV",
    "ETAP_02_status_741_bialek_panelu.csv"
  ),
  show_col_types = FALSE,
  locale = readr::locale(
    decimal_mark = ",",
    grouping_mark = " ",
    encoding = "UTF-8"
  )
) %>%
  mutate(Bialko = stringr::str_trim(as.character(Bialko))) %>%
  filter(!is.na(Bialko), nzchar(Bialko)) %>%
  distinct(Bialko, .keep_all = TRUE)

panel_ids <- status_panel$Bialko

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

podsumowanie_gradienty_etap02 <- readr::read_csv2(
  file.path(folder_etap02, "04_CSV", "ETAP_02_porownanie_0_vs_1_DPMC_30_60_szczegolowe.csv"),
  show_col_types = FALSE,
  locale = readr::locale(decimal_mark = ",", grouping_mark = " ", encoding = "UTF-8")
) %>%
  mutate(
    Etap_ID = as.integer(na_liczbe(Etap_ID)),
    Liczba_bialek = as.integer(na_liczbe(Liczba_bialek)),
    Liczba_peptydow = as.integer(na_liczbe(Liczba_peptydow))
  )

oczekiwane_top <- function(limit_id, kolumna) {
  x <- podsumowanie_top_etap02 %>%
    filter(Zbior == "Panel 741 białek", Limit_ID == limit_id) %>%
    pull(all_of(kolumna))
  if (length(x) != 1L) stop("Nie udało się odczytać liczebności Top2/Top3 z ETAPU 02.")
  x[[1]]
}

oczekiwane_gradient <- function(gradient, etap_id, kolumna) {
  x <- podsumowanie_gradienty_etap02 %>%
    filter(
      Wariant_ID == "1DPMC",
      Gradient == gradient,
      Zbior == "Panel 741 białek",
      Etap_ID == etap_id
    ) %>%
    pull(all_of(kolumna))
  if (length(x) != 1L) stop("Nie udało się odczytać liczebności gradientu z ETAPU 02.")
  x[[1]]
}

# Wczytuje bazową listę Top 2/Top 3 z etapu 03 i ujednolica jej strukturę.
wczytaj_liste_ETAP03 <- function(nazwa_pliku, typ_listy, maks_n) {
  df <- readr::read_csv2(
    file.path(
      folder_etap03,
      "02_CSV_do_dalszego_przygotowania",
      nazwa_pliku
    ),
    show_col_types = FALSE,
    locale = readr::locale(
      decimal_mark = ",",
      grouping_mark = " ",
      encoding = "UTF-8"
    ),
    na = c("", "NA", "NaN")
  )

  wymagane <- c(
    "Bialko",
    "Nazwa_bialka",
    "Peptyd",
    "Wybrany_ladunek",
    "N_plikow_zmierzony",
    "RT_sredni",
    "RT_start_sredni",
    "RT_stop_sredni",
    "Najlepszy_PEP",
    "Max_intensywnosc",
    "Wybrane_mz",
    "Max_dobrych_fragmentow",
    "Max_srednia_top_fragmentow",
    "Max_suma_top_fragmentow",
    "Najmocniejszy_fragment",
    "Max_sredni_score_top3",
    "Max_min_score_top3",
    "Max_udzial_runow_top3",
    "Pozycja_w_bialku"
  )

  sprawdz_kolumny(df, wymagane, typ_listy)

  df %>%
    select(-matches("QValue|Q[.]Value|Q_Value", ignore.case = TRUE)) %>%
    mutate(
      Bialko = as.character(Bialko),
      Nazwa_bialka = as.character(Nazwa_bialka),
      Peptyd = as.character(Peptyd),
      Wybrany_ladunek = as.integer(na_liczbe(Wybrany_ladunek)),
      N_plikow_zmierzony = as.integer(na_liczbe(N_plikow_zmierzony)),
      RT_sredni = na_liczbe(RT_sredni),
      RT_start_sredni = na_liczbe(RT_start_sredni),
      RT_stop_sredni = na_liczbe(RT_stop_sredni),
      Najlepszy_PEP = na_liczbe(Najlepszy_PEP),
      Max_intensywnosc = na_liczbe(Max_intensywnosc),
      Wybrane_mz = na_liczbe(Wybrane_mz),
      Max_dobrych_fragmentow = as.integer(
        na_liczbe(Max_dobrych_fragmentow)
      ),
      Max_srednia_top_fragmentow = na_liczbe(
        Max_srednia_top_fragmentow
      ),
      Max_suma_top_fragmentow = na_liczbe(
        Max_suma_top_fragmentow
      ),
      Najmocniejszy_fragment = na_liczbe(Najmocniejszy_fragment),
      Max_sredni_score_top3 = na_liczbe(Max_sredni_score_top3),
      Max_min_score_top3 = na_liczbe(Max_min_score_top3),
      Max_udzial_runow_top3 = na_liczbe(Max_udzial_runow_top3),
      Pozycja_w_bialku = as.integer(na_liczbe(Pozycja_w_bialku)),
      Typ_listy = typ_listy,
      Maks_peptydow_na_bialko = maks_n
    ) %>%
    distinct(Bialko, Peptyd, .keep_all = TRUE) %>%
    arrange(Bialko, Pozycja_w_bialku)
}

lista_top2 <- wczytaj_liste_ETAP03(
  "CSV_02A_PANEL_741_lista_peptydowo_prekursorowa_top2_PO_RAGGED.csv",
  "Top 2",
  2L
)

lista_top3 <- wczytaj_liste_ETAP03(
  "CSV_02B_PANEL_741_lista_peptydowo_prekursorowa_top3_PO_RAGGED.csv",
  "Top 3",
  3L
)

obserwowane_30_panel <- obserwowane_30_all %>%
  filter(Bialko %in% panel_ids)

obserwowane_60_panel <- obserwowane_60_all %>%
  filter(Bialko %in% panel_ids)

kwalifikowane_30_panel <- kwalifikowane_30_all %>%
  filter(Bialko %in% panel_ids)

kwalifikowane_60_panel <- kwalifikowane_60_all %>%
  filter(Bialko %in% panel_ids)

jakosc_30_panel <- jakosc_30_all %>%
  filter(Bialko %in% panel_ids)

jakosc_60_panel <- jakosc_60_all %>%
  filter(Bialko %in% panel_ids)


# ------------------------------------------------------------------------------
# 7. WALIDACJA PRZENOSZENIA RT60 -> RT30
# ------------------------------------------------------------------------------

# Walidacja jest wykonywana tylko dla peptydów obserwowanych w obu gradientach.

wspolne_RT <- jakosc_30_all %>%
  transmute(
    Bialko,
    Nazwa_bialka,
    Peptyd,
    RT30_obserwowany = RT_sredni,
    PEP_30 = Najlepszy_PEP,
    Fragmenty_30 = Max_dobrych_fragmentow
  ) %>%
  inner_join(
    jakosc_60_all %>%
      transmute(
        Bialko,
        Peptyd,
        RT60_obserwowany = RT_sredni,
        PEP_60 = Najlepszy_PEP,
        Fragmenty_60 = Max_dobrych_fragmentow
      ),
    by = c("Bialko", "Peptyd")
  ) %>%
  mutate(
    iRT_z_RT60 = RT60_na_iRT(RT60_obserwowany),
    RT30_przewidziany = iRT_na_RT30(iRT_z_RT60),
    Blad_predykcji_min = RT30_przewidziany - RT30_obserwowany,
    Blad_bezwzgledny_min = abs(Blad_predykcji_min),
    W_panelu_741 = Bialko %in% panel_ids
  ) %>%
  filter(
    is.finite(RT30_obserwowany),
    is.finite(RT60_obserwowany),
    is.finite(RT30_przewidziany),
    is.finite(Blad_predykcji_min)
  )

# Oblicza miary błędu transferu czasu retencji pomiędzy gradientami.
podsumuj_blad <- function(df, zakres) {
  if (nrow(df) == 0) {
    return(tibble(
      Zakres = zakres,
      N_peptydow = 0L,
      Blad_sredni_min = NA_real_,
      Blad_mediana_min = NA_real_,
      MAE_min = NA_real_,
      RMSE_min = NA_real_,
      Blad_abs_90_percentyl_min = NA_real_,
      Blad_abs_95_percentyl_min = NA_real_,
      R2_obserwowany_vs_przewidywany = NA_real_
    ))
  }

  tibble(
    Zakres = zakres,
    N_peptydow = nrow(df),
    Blad_sredni_min = mean(df$Blad_predykcji_min, na.rm = TRUE),
    Blad_mediana_min = median(df$Blad_predykcji_min, na.rm = TRUE),
    MAE_min = mean(df$Blad_bezwzgledny_min, na.rm = TRUE),
    RMSE_min = sqrt(mean(df$Blad_predykcji_min^2, na.rm = TRUE)),
    Blad_abs_90_percentyl_min = as.numeric(
      quantile(df$Blad_bezwzgledny_min, 0.90, na.rm = TRUE)
    ),
    Blad_abs_95_percentyl_min = as.numeric(
      quantile(df$Blad_bezwzgledny_min, 0.95, na.rm = TRUE)
    ),
    R2_obserwowany_vs_przewidywany = if (nrow(df) >= 2) {
      cor(
        df$RT30_obserwowany,
        df$RT30_przewidziany,
        use = "complete.obs"
      )^2
    } else {
      NA_real_
    }
  )
}

podsumowanie_bledu <- bind_rows(
  podsumuj_blad(
    wspolne_RT,
    "Cała FASTA human-only po PEP i fragmentach"
  ),
  podsumuj_blad(
    wspolne_RT %>% filter(W_panelu_741),
    "Panel 741 po PEP i fragmentach"
  )
)

margines_95_panel <- podsumowanie_bledu %>%
  filter(Zakres == "Panel 741 po PEP i fragmentach") %>%
  pull(Blad_abs_95_percentyl_min)

margines_95_all <- podsumowanie_bledu %>%
  filter(Zakres == "Cała FASTA human-only po PEP i fragmentach") %>%
  pull(Blad_abs_95_percentyl_min)

margines_95_uzyty <- if (
  length(margines_95_panel) == 1 &&
    is.finite(margines_95_panel) &&
    sum(wspolne_RT$W_panelu_741) >= 100
) {
  margines_95_panel
} else {
  margines_95_all
}


# ------------------------------------------------------------------------------
# 8. BUDOWA PUL I LIST ROZSZERZONYCH
# ------------------------------------------------------------------------------

# Najpierw wykorzystywane są kandydaty z obserwowanym RT30, a dopiero później peptydy wymagające przewidywania RT30.

rt30_lookup <- obserwowane_30_panel %>%
  filter(
    is.finite(RT_sredni),
    is.finite(RT_start_sredni),
    is.finite(RT_stop_sredni)
  ) %>%
  transmute(
    Bialko,
    Peptyd,
    RT30_obserwowany = RT_sredni,
    RT30_start_obserwowany = RT_start_sredni,
    RT30_stop_obserwowany = RT_stop_sredni
  ) %>%
  distinct(Bialko, Peptyd, .keep_all = TRUE)

rt60_lookup <- obserwowane_60_panel %>%
  filter(is.finite(RT_sredni)) %>%
  transmute(
    Bialko,
    Peptyd,
    RT60_obserwowany_lookup = RT_sredni
  ) %>%
  distinct(Bialko, Peptyd, .keep_all = TRUE)

# Uzupełnia wolne miejsca na liście bazowej zgodnie z kolejnością źródeł RT i regułami jakości.
zbuduj_wariant <- function(lista_bazowa, maks_n, typ_listy) {
  klucze_bazy <- lista_bazowa %>%
    select(Bialko, Peptyd)

  # Peptydy zaobserwowane w 30 min, lecz jakość pobierana z 60 min.
  pula_RT30_obserwowany_jakosc60 <- jakosc_60_panel %>%
    anti_join(klucze_bazy, by = c("Bialko", "Peptyd")) %>%
    inner_join(rt30_lookup, by = c("Bialko", "Peptyd")) %>%
    ranking_kandydatow()

  # Prawdziwe „tylko 60 min”: anti-join względem WSZYSTKICH peptydów
  # zaobserwowanych w 30 min, a nie tylko względem puli po PEP/fragmentach.
  pula_prawdziwie_tylko_60 <- jakosc_60_panel %>%
    anti_join(
      obserwowane_30_panel %>% select(Bialko, Peptyd),
      by = c("Bialko", "Peptyd")
    ) %>%
    anti_join(klucze_bazy, by = c("Bialko", "Peptyd")) %>%
    ranking_kandydatow()

  liczby_bazowe <- tibble(Bialko = panel_ids) %>%
    left_join(
      lista_bazowa %>%
        count(Bialko, name = "Liczba_bazowa"),
      by = "Bialko"
    ) %>%
    mutate(
      Liczba_bazowa = replace_na(Liczba_bazowa, 0L),
      Wolne_miejsca = pmax(maks_n - Liczba_bazowa, 0L)
    )

  wybrane_RT30_obserwowany <- pula_RT30_obserwowany_jakosc60 %>%
    inner_join(liczby_bazowe, by = "Bialko") %>%
    filter(
      Wolne_miejsca > 0,
      Ranking_w_puli <= Wolne_miejsca
    ) %>%
    mutate(
      Kolejnosc_dodania = Ranking_w_puli,
      Typ_dodania = "RT30_obserwowany_jakosc_60min"
    )

  liczby_dodane_obserwowane <- wybrane_RT30_obserwowany %>%
    count(Bialko, name = "Dodane_RT30_obserwowany")

  pozostale_miejsca <- liczby_bazowe %>%
    left_join(liczby_dodane_obserwowane, by = "Bialko") %>%
    mutate(
      Dodane_RT30_obserwowany = replace_na(
        Dodane_RT30_obserwowany,
        0L
      ),
      Pozostale_miejsca = pmax(
        Wolne_miejsca - Dodane_RT30_obserwowany,
        0L
      )
    )

  wybrane_prawdziwie_tylko_60 <- pula_prawdziwie_tylko_60 %>%
    inner_join(pozostale_miejsca, by = "Bialko") %>%
    filter(
      Pozostale_miejsca > 0,
      Ranking_w_puli <= Pozostale_miejsca
    ) %>%
    mutate(
      Kolejnosc_dodania = Dodane_RT30_obserwowany + Ranking_w_puli,
      Typ_dodania = "prawdziwie_tylko_60min_predykcja_iRT"
    )

  baza_std <- lista_bazowa %>%
    left_join(rt60_lookup, by = c("Bialko", "Peptyd")) %>%
    transmute(
      Bialko,
      Nazwa_bialka,
      Peptyd,
      Wybrany_ladunek,
      Wybrane_mz,
      N_plikow_zmierzony,
      Najlepszy_PEP,
      Max_dobrych_fragmentow,
      Max_srednia_top_fragmentow,
      Max_suma_top_fragmentow,
      Max_intensywnosc,
      Najmocniejszy_fragment,
      Max_sredni_score_top3,
      Max_min_score_top3,
      Max_udzial_runow_top3,
      Pozycja_w_bialku = Pozycja_w_bialku,
      Typ_listy = typ_listy,
      Maks_peptydow_na_bialko = maks_n,
      Typ_rekordu = "bazowy_ETAP03",
      Zrodlo_RT = "DIA-NN_30min_obserwowany",
      Zrodlo_jakosci = "DIA-NN_30min_EXPORT-QUANT",
      Zrodlo_ladunku_i_mz = "ladunek_DIA-NN_30min; m/z_teoretyczne_CAM_C",
      Czy_obserwowany_w_30min = TRUE,
      Czy_prawdziwie_tylko_60min = FALSE,
      Czy_RT_przewidywany = FALSE,
      RT30_obserwowany = RT_sredni,
      RT30_start_obserwowany = RT_start_sredni,
      RT30_stop_obserwowany = RT_stop_sredni,
      RT60_obserwowany = RT60_obserwowany_lookup,
      iRT = RT30_na_iRT(RT_sredni),
      RT30_przewidywany = NA_real_,
      RT30_docelowy = RT_sredni,
      RT30_zakres_orientacyjny_min = RT_start_sredni,
      RT30_zakres_orientacyjny_max = RT_stop_sredni,
      Margines_bledu_95_min = NA_real_,
      Uwagi = paste(
        "Peptyd bazowy z ETAPU 03.",
        "RT i parametry jakości pochodzą z gradientu 30-minutowego."
      ),
      ID_peptydu = paste(Bialko, Peptyd, sep = "|"),
      ID_prekursora = paste(Bialko, Peptyd, Wybrany_ladunek, sep = "|")
    )

  dodane_obserwowane_std <- wybrane_RT30_obserwowany %>%
    transmute(
      Bialko,
      Nazwa_bialka,
      Peptyd,
      Wybrany_ladunek,
      Wybrane_mz,
      N_plikow_zmierzony,
      Najlepszy_PEP,
      Max_dobrych_fragmentow,
      Max_srednia_top_fragmentow,
      Max_suma_top_fragmentow,
      Max_intensywnosc,
      Najmocniejszy_fragment,
      Max_sredni_score_top3,
      Max_min_score_top3,
      Max_udzial_runow_top3,
      Pozycja_w_bialku = Liczba_bazowa + Kolejnosc_dodania,
      Typ_listy = typ_listy,
      Maks_peptydow_na_bialko = maks_n,
      Typ_rekordu = "dodany_RT30_obserwowany_jakosc60",
      Zrodlo_RT = "DIA-NN_30min_obserwowany",
      Zrodlo_jakosci = "DIA-NN_60min_EXPORT-QUANT",
      Zrodlo_ladunku_i_mz = "ladunek_DIA-NN_60min; m/z_teoretyczne_CAM_C",
      Czy_obserwowany_w_30min = TRUE,
      Czy_prawdziwie_tylko_60min = FALSE,
      Czy_RT_przewidywany = FALSE,
      RT30_obserwowany,
      RT30_start_obserwowany,
      RT30_stop_obserwowany,
      RT60_obserwowany = RT_sredni,
      iRT = RT30_na_iRT(RT30_obserwowany),
      RT30_przewidywany = NA_real_,
      RT30_docelowy = RT30_obserwowany,
      RT30_zakres_orientacyjny_min = RT30_start_obserwowany,
      RT30_zakres_orientacyjny_max = RT30_stop_obserwowany,
      Margines_bledu_95_min = NA_real_,
      Uwagi = paste(
        "Peptyd był zaobserwowany w 30 min, dlatego użyto rzeczywistego RT30.",
        "PEP, kwalifikujące się fragmenty, Fr.Score, intensywności, ładunek i m/z pochodzą z 60 min."
      ),
      ID_peptydu = paste(Bialko, Peptyd, sep = "|"),
      ID_prekursora = paste(Bialko, Peptyd, Wybrany_ladunek, sep = "|")
    )

  dodane_predykcja_std <- wybrane_prawdziwie_tylko_60 %>%
    mutate(
      iRT_obliczony = RT60_na_iRT(RT_sredni),
      RT30_pred = iRT_na_RT30(iRT_obliczony)
    ) %>%
    transmute(
      Bialko,
      Nazwa_bialka,
      Peptyd,
      Wybrany_ladunek,
      Wybrane_mz,
      N_plikow_zmierzony,
      Najlepszy_PEP,
      Max_dobrych_fragmentow,
      Max_srednia_top_fragmentow,
      Max_suma_top_fragmentow,
      Max_intensywnosc,
      Najmocniejszy_fragment,
      Max_sredni_score_top3,
      Max_min_score_top3,
      Max_udzial_runow_top3,
      Pozycja_w_bialku = Liczba_bazowa + Kolejnosc_dodania,
      Typ_listy = typ_listy,
      Maks_peptydow_na_bialko = maks_n,
      Typ_rekordu = "dodany_prawdziwie_tylko60_predykcja_iRT",
      Zrodlo_RT = "przewidywany_z_RT60_przez_iRT",
      Zrodlo_jakosci = "DIA-NN_60min_EXPORT-QUANT",
      Zrodlo_ladunku_i_mz = "ladunek_DIA-NN_60min; m/z_teoretyczne_CAM_C",
      Czy_obserwowany_w_30min = FALSE,
      Czy_prawdziwie_tylko_60min = TRUE,
      Czy_RT_przewidywany = TRUE,
      RT30_obserwowany = NA_real_,
      RT30_start_obserwowany = NA_real_,
      RT30_stop_obserwowany = NA_real_,
      RT60_obserwowany = RT_sredni,
      iRT = iRT_obliczony,
      RT30_przewidywany = RT30_pred,
      RT30_docelowy = RT30_pred,
      RT30_zakres_orientacyjny_min = RT30_pred - margines_95_uzyty,
      RT30_zakres_orientacyjny_max = RT30_pred + margines_95_uzyty,
      Margines_bledu_95_min = margines_95_uzyty,
      Uwagi = paste(
        "Peptyd nie był zaobserwowany w żadnym pliku 30-minutowym.",
        "RT30 przewidziano z obserwowanego RT60 przez iRT.",
        "Zakres +/- P95 jest orientacyjny i nie stanowi finalnego okna scheduled PRM."
      ),
      ID_peptydu = paste(Bialko, Peptyd, sep = "|"),
      ID_prekursora = paste(Bialko, Peptyd, Wybrany_ladunek, sep = "|")
    )

  lista_rozszerzona <- bind_rows(
    baza_std,
    dodane_obserwowane_std,
    dodane_predykcja_std
  ) %>%
    arrange(Bialko, Pozycja_w_bialku, Najlepszy_PEP, Peptyd)

  counts_baza <- baza_std %>%
    count(Bialko, name = "Peptydy_bazowe")

  counts_obs <- dodane_obserwowane_std %>%
    count(Bialko, name = "Peptydy_dodane_RT30_obserwowany")

  counts_pred <- dodane_predykcja_std %>%
    count(Bialko, name = "Peptydy_dodane_RT30_przewidywany")

  status_wariant <- status_panel %>%
    select(Bialko, everything()) %>%
    left_join(counts_baza, by = "Bialko") %>%
    left_join(counts_obs, by = "Bialko") %>%
    left_join(counts_pred, by = "Bialko") %>%
    mutate(
      Peptydy_bazowe = replace_na(Peptydy_bazowe, 0L),
      Peptydy_dodane_RT30_obserwowany = replace_na(
        Peptydy_dodane_RT30_obserwowany,
        0L
      ),
      Peptydy_dodane_RT30_przewidywany = replace_na(
        Peptydy_dodane_RT30_przewidywany,
        0L
      ),
      Peptydy_razem = Peptydy_bazowe +
        Peptydy_dodane_RT30_obserwowany +
        Peptydy_dodane_RT30_przewidywany,
      Typ_listy = typ_listy,
      Maks_peptydow_na_bialko = maks_n,
      Status_ETAP05 = case_when(
        Peptydy_bazowe == maks_n ~
          paste0("Lista bazowa pełna (", maks_n, ")"),
        Peptydy_bazowe > 0 &
          Peptydy_dodane_RT30_obserwowany > 0 &
          Peptydy_dodane_RT30_przewidywany > 0 ~
          "Uzupełniono RT30 obserwowanym i predykcją",
        Peptydy_bazowe > 0 &
          Peptydy_dodane_RT30_obserwowany > 0 ~
          "Uzupełniono peptydem z RT30 obserwowanym; jakość 60 min",
        Peptydy_bazowe > 0 &
          Peptydy_dodane_RT30_przewidywany > 0 ~
          "Uzupełniono prawdziwym tylko-60; RT30 przewidywany",
        Peptydy_bazowe == 0 &
          Peptydy_dodane_RT30_obserwowany > 0 &
          Peptydy_dodane_RT30_przewidywany > 0 ~
          "Białko dodane przez RT30 obserwowany i predykcję",
        Peptydy_bazowe == 0 &
          Peptydy_dodane_RT30_obserwowany > 0 ~
          "Białko dodane z RT30 obserwowanym; jakość 60 min",
        Peptydy_bazowe == 0 &
          Peptydy_dodane_RT30_przewidywany > 0 ~
          "Białko dodane prawdziwym tylko-60; RT30 przewidywany",
        Peptydy_bazowe > 0 ~
          paste0(
            "Pozostało ",
            Peptydy_bazowe,
            "/",
            maks_n,
            "; brak dalszego kandydata"
          ),
        TRUE ~
          "Brak peptydu na liście rozszerzonej"
      )
    )

  list(
    typ = typ_listy,
    maks_n = maks_n,
    baza = baza_std,
    pula_obserwowany = pula_RT30_obserwowany_jakosc60,
    pula_tylko60 = pula_prawdziwie_tylko_60,
    dodane_obserwowany = dodane_obserwowane_std,
    dodane_predykcja = dodane_predykcja_std,
    rozszerzona = lista_rozszerzona,
    status = status_wariant
  )
}

wariant_top2 <- zbuduj_wariant(lista_top2, 2L, "Top 2")
wariant_top3 <- zbuduj_wariant(lista_top3, 3L, "Top 3")


# ------------------------------------------------------------------------------
# 9. PODSUMOWANIA
# ------------------------------------------------------------------------------

podsumuj_wariant <- function(w) {
  bialka_baza <- unique(w$baza$Bialko)
  bialka_ext <- unique(w$rozszerzona$Bialko)

  bind_rows(
    tibble(
      Typ_listy = w$typ,
      Wariant = "Bazowa ETAP 03",
      Liczba_bialek = n_distinct(w$baza$Bialko),
      Liczba_peptydow = nrow(w$baza),
      RT30_obserwowany = nrow(w$baza),
      RT30_przewidywany = 0L,
      Jakosc_z_30min = nrow(w$baza),
      Jakosc_z_60min = 0L,
      Bialka_dodane_wzgledem_bazy = 0L
    ),
    tibble(
      Typ_listy = w$typ,
      Wariant = "Rozszerzona ETAP 05",
      Liczba_bialek = n_distinct(w$rozszerzona$Bialko),
      Liczba_peptydow = nrow(w$rozszerzona),
      RT30_obserwowany = sum(!w$rozszerzona$Czy_RT_przewidywany),
      RT30_przewidywany = sum(w$rozszerzona$Czy_RT_przewidywany),
      Jakosc_z_30min = sum(w$rozszerzona$Zrodlo_jakosci == "DIA-NN_30min_EXPORT-QUANT"),
      Jakosc_z_60min = sum(w$rozszerzona$Zrodlo_jakosci == "DIA-NN_60min_EXPORT-QUANT"),
      Bialka_dodane_wzgledem_bazy = length(setdiff(bialka_ext, bialka_baza))
    )
  )
}

podsumowanie_list <- bind_rows(
  podsumuj_wariant(wariant_top2),
  podsumuj_wariant(wariant_top3)
)

dystrybucja_peptydow <- bind_rows(
  wariant_top2$rozszerzona %>%
    count(Bialko, name = "Liczba_peptydow_na_bialko") %>%
    count(Liczba_peptydow_na_bialko, name = "Liczba_bialek") %>%
    mutate(Typ_listy = "Top 2"),
  wariant_top3$rozszerzona %>%
    count(Bialko, name = "Liczba_peptydow_na_bialko") %>%
    count(Liczba_peptydow_na_bialko, name = "Liczba_bialek") %>%
    mutate(Typ_listy = "Top 3")
) %>%
  group_by(Typ_listy) %>%
  mutate(
    Procent_bialek = round(100 * Liczba_bialek / sum(Liczba_bialek), 2)
  ) %>%
  ungroup()

podsumowanie_statusow <- bind_rows(
  wariant_top2$status,
  wariant_top3$status
) %>%
  count(Typ_listy, Status_ETAP05, name = "Liczba_bialek") %>%
  group_by(Typ_listy) %>%
  mutate(
    Procent_panelu = round(100 * Liczba_bialek / sum(Liczba_bialek), 2)
  ) %>%
  ungroup() %>%
  arrange(Typ_listy, desc(Liczba_bialek))


# ------------------------------------------------------------------------------
# 10. KONTROLE POPRAWNOŚCI
# ------------------------------------------------------------------------------

# Kontrole sprawdzają między innymi źródło RT30, liczebności list, limity peptydów na białko i zachowanie kryteriów jakości.

klucz <- c("Bialko", "Peptyd")

top2_nie_w_top3 <- wariant_top2$rozszerzona %>%
  anti_join(
    wariant_top3$rozszerzona %>% select(all_of(klucz)),
    by = klucz
  )

predykcja_obecna_w_30_top2 <- wariant_top2$dodane_predykcja %>%
  inner_join(
    obserwowane_30_panel %>% select(all_of(klucz)),
    by = klucz
  )

predykcja_obecna_w_30_top3 <- wariant_top3$dodane_predykcja %>%
  inner_join(
    obserwowane_30_panel %>% select(all_of(klucz)),
    by = klucz
  )

# Wykonuje zestaw kontroli dla gotowej listy rozszerzonej.
sprawdz_wariant <- function(w) {
  duplikaty_peptydow <- w$rozszerzona %>%
    count(Bialko, Peptyd) %>%
    filter(n > 1)

  duplikaty_prekursorow <- w$rozszerzona %>%
    count(ID_prekursora) %>%
    filter(n > 1)

  nieunikatowe <- w$rozszerzona %>%
    distinct(Bialko, Peptyd) %>%
    count(Peptyd, name = "Liczba_bialek") %>%
    filter(Liczba_bialek > 1)

  liczby_na_bialko <- w$rozszerzona %>%
    count(Bialko, name = "n")

  baza_niezachowana <- w$baza %>%
    anti_join(
      w$rozszerzona %>% select(Bialko, Peptyd, Typ_rekordu),
      by = c("Bialko", "Peptyd")
    )

  dodane_obserwowane_juz_w_jakosci30 <- w$dodane_obserwowany %>%
    inner_join(
      jakosc_30_panel %>% select(Bialko, Peptyd),
      by = c("Bialko", "Peptyd")
    )

  bind_rows(
    kontrola(
      paste0(w$typ, ": brak duplikatów białko-peptyd"),
      nrow(duplikaty_peptydow) == 0,
      nrow(duplikaty_peptydow),
      0
    ),
    kontrola(
      paste0(w$typ, ": brak duplikatów prekursorów"),
      nrow(duplikaty_prekursorow) == 0,
      nrow(duplikaty_prekursorow),
      0
    ),
    kontrola(
      paste0(w$typ, ": globalna unikatowość sekwencji"),
      nrow(nieunikatowe) == 0,
      nrow(nieunikatowe),
      0
    ),
    kontrola(
      paste0(w$typ, ": maksymalna liczba peptydów na białko"),
      max(liczby_na_bialko$n) <= w$maks_n,
      max(liczby_na_bialko$n),
      paste0("<= ", w$maks_n)
    ),
    kontrola(
      paste0(w$typ, ": zachowano wszystkie rekordy bazowe"),
      nrow(baza_niezachowana) == 0,
      nrow(baza_niezachowana),
      0
    ),
    kontrola(
      paste0(w$typ, ": PEP <= 0,01"),
      all(
        is.finite(w$rozszerzona$Najlepszy_PEP) &
          w$rozszerzona$Najlepszy_PEP <= prog_PEP
      ),
      sum(
        !is.finite(w$rozszerzona$Najlepszy_PEP) |
          w$rozszerzona$Najlepszy_PEP > prog_PEP
      ),
      0
    ),
    kontrola(
      paste0(w$typ, ": co najmniej 3 kwalifikujące się fragmenty"),
      all(
        !is.na(w$rozszerzona$Max_dobrych_fragmentow) &
          w$rozszerzona$Max_dobrych_fragmentow >= min_dobrych_fragmentow
      ),
      sum(
        is.na(w$rozszerzona$Max_dobrych_fragmentow) |
          w$rozszerzona$Max_dobrych_fragmentow < min_dobrych_fragmentow
      ),
      0
    ),
    kontrola(
      paste0(w$typ, ": metryki Fr.Score top3 są skończone"),
      all(
        is.finite(w$rozszerzona$Max_sredni_score_top3) &
          is.finite(w$rozszerzona$Max_min_score_top3) &
          is.finite(w$rozszerzona$Max_udzial_runow_top3)
      ),
      sum(
        !is.finite(w$rozszerzona$Max_sredni_score_top3) |
          !is.finite(w$rozszerzona$Max_min_score_top3) |
          !is.finite(w$rozszerzona$Max_udzial_runow_top3)
      ),
      0
    ),
    kontrola(
      paste0(w$typ, ": wyłącznie ładunki 2+ i 3+"),
      all(w$rozszerzona$Wybrany_ladunek %in% dozwolone_ladunki),
      sum(!w$rozszerzona$Wybrany_ladunek %in% dozwolone_ladunki),
      0
    ),
    kontrola(
      paste0(w$typ, ": brak ragged"),
      !any(czy_ragged(w$rozszerzona$Peptyd)),
      sum(czy_ragged(w$rozszerzona$Peptyd)),
      0
    ),
    kontrola(
      paste0(w$typ, ": brak metioniny"),
      !any(stringr::str_detect(w$rozszerzona$Peptyd, "M")),
      sum(stringr::str_detect(w$rozszerzona$Peptyd, "M")),
      0
    ),
    kontrola(
      paste0(w$typ, ": brak standardu iRT"),
      !any(czy_standard_iRT(w$rozszerzona$Bialko)),
      sum(czy_standard_iRT(w$rozszerzona$Bialko)),
      0
    ),
    kontrola(
      paste0(w$typ, ": tylko białka panelu 741"),
      all(w$rozszerzona$Bialko %in% panel_ids),
      sum(!w$rozszerzona$Bialko %in% panel_ids),
      0
    ),
    kontrola(
      paste0(w$typ, ": RT przewidywany tylko dla prawdziwie tylko-60"),
      all(
        w$rozszerzona$Czy_RT_przewidywany ==
          w$rozszerzona$Czy_prawdziwie_tylko_60min
      ),
      sum(
        w$rozszerzona$Czy_RT_przewidywany !=
          w$rozszerzona$Czy_prawdziwie_tylko_60min
      ),
      0
    ),
    kontrola(
      paste0(w$typ, ": obserwowany RT30 nie został zastąpiony predykcją"),
      all(
        !w$rozszerzona$Czy_obserwowany_w_30min |
          !w$rozszerzona$Czy_RT_przewidywany
      ),
      sum(
        w$rozszerzona$Czy_obserwowany_w_30min &
          w$rozszerzona$Czy_RT_przewidywany
      ),
      0
    ),
    kontrola(
      paste0(w$typ, ": dla obserwowanych RT30 docelowy jest równy obserwowanemu"),
      all(
        !w$rozszerzona$Czy_obserwowany_w_30min |
          abs(
            w$rozszerzona$RT30_docelowy -
              w$rozszerzona$RT30_obserwowany
          ) < 1e-10
      ),
      sum(
        w$rozszerzona$Czy_obserwowany_w_30min &
          abs(
            w$rozszerzona$RT30_docelowy -
              w$rozszerzona$RT30_obserwowany
          ) >= 1e-10
      ),
      0
    ),
    kontrola(
      paste0(w$typ, ": dodatki z jakością 60 min nie należą do puli jakości 30 min"),
      nrow(dodane_obserwowane_juz_w_jakosci30) == 0,
      nrow(dodane_obserwowane_juz_w_jakosci30),
      0
    ),
    kontrola(
      paste0(w$typ, ": przewidywane rekordy nie mają obserwowanego RT30"),
      all(
        !w$rozszerzona$Czy_RT_przewidywany |
          is.na(w$rozszerzona$RT30_obserwowany)
      ),
      sum(
        w$rozszerzona$Czy_RT_przewidywany &
          !is.na(w$rozszerzona$RT30_obserwowany)
      ),
      0
    ),
    kontrola(
      paste0(w$typ, ": wszystkie docelowe RT30 są skończone"),
      all(is.finite(w$rozszerzona$RT30_docelowy)),
      sum(!is.finite(w$rozszerzona$RT30_docelowy)),
      0
    ),
    kontrola(
      paste0(w$typ, ": dokładnie 741 statusów"),
      nrow(w$status) == length(panel_ids) &&
        n_distinct(w$status$Bialko) == length(panel_ids),
      paste(nrow(w$status), n_distinct(w$status$Bialko), sep = "; "),
      paste(length(panel_ids), length(panel_ids), sep = "; ")
    ),
    kontrola(
      paste0(w$typ, ": brak Q.Value w liście końcowej"),
      !any(
        stringr::str_detect(
          colnames(w$rozszerzona),
          stringr::regex("QValue|Q[.]Value|Q_Value", ignore_case = TRUE)
        )
      ),
      paste(
        colnames(w$rozszerzona)[
          stringr::str_detect(
            colnames(w$rozszerzona),
            stringr::regex("QValue|Q[.]Value|Q_Value", ignore_case = TRUE)
          )
        ],
        collapse = "; "
      ),
      "brak"
    )
  )
}

kontrole <- bind_rows(
  kontrola(
    "Kalibracja zawiera gradient 30min",
    "30min" %in% kalibracja$Gradient,
    "30min" %in% kalibracja$Gradient,
    TRUE
  ),
  kontrola(
    "Kalibracja zawiera gradient 60min",
    "60min" %in% kalibracja$Gradient,
    "60min" %in% kalibracja$Gradient,
    TRUE
  ),
  kontrola(
    "Współczynniki kalibracji są skończone i dodatnie",
    all(is.finite(c(a30, b30, a60, b60))) && a30 > 0 && a60 > 0,
    paste(a30, b30, a60, b60, sep = "; "),
    "skończone; a30 i a60 > 0"
  ),
  kontrola(
    "R2 kalibracji >= 0,90",
    all(kalibracja$R2 >= 0.90),
    paste(round(kalibracja$R2, 6), collapse = "; "),
    ">= 0,90"
  ),
  kontrola(
    "Panel zawiera 741 unikatowych białek",
    length(panel_ids) == 741L,
    length(panel_ids),
    741
  ),
  kontrola(
    "Top 2 bazowa — liczba peptydów zgodna z ETAPEM 02",
    nrow(wariant_top2$baza) == oczekiwane_top(2L, "Liczba_peptydow"),
    nrow(wariant_top2$baza),
    oczekiwane_top(2L, "Liczba_peptydow")
  ),
  kontrola(
    "Top 2 bazowa — liczba białek zgodna z ETAPEM 02",
    n_distinct(wariant_top2$baza$Bialko) == oczekiwane_top(2L, "Liczba_bialek"),
    n_distinct(wariant_top2$baza$Bialko),
    oczekiwane_top(2L, "Liczba_bialek")
  ),
  kontrola(
    "Top 3 bazowa — liczba peptydów zgodna z ETAPEM 02",
    nrow(wariant_top3$baza) == oczekiwane_top(3L, "Liczba_peptydow"),
    nrow(wariant_top3$baza),
    oczekiwane_top(3L, "Liczba_peptydow")
  ),
  kontrola(
    "Top 3 bazowa — liczba białek zgodna z ETAPEM 02",
    n_distinct(wariant_top3$baza$Bialko) == oczekiwane_top(3L, "Liczba_bialek"),
    n_distinct(wariant_top3$baza$Bialko),
    oczekiwane_top(3L, "Liczba_bialek")
  ),
  kontrola(
    "Panel 30 min — liczba rekordów z RT/ładunkiem zgodna z ETAPEM 02",
    nrow(kwalifikowane_30_panel) == oczekiwane_gradient("30 min", 3L, "Liczba_peptydow"),
    nrow(kwalifikowane_30_panel),
    oczekiwane_gradient("30 min", 3L, "Liczba_peptydow")
  ),
  kontrola(
    "Panel 30 min — jakość PEP + fragmenty zgodna z ETAPEM 02",
    nrow(jakosc_30_panel) == oczekiwane_gradient("30 min", 5L, "Liczba_peptydow"),
    nrow(jakosc_30_panel),
    oczekiwane_gradient("30 min", 5L, "Liczba_peptydow")
  ),
  kontrola(
    "Panel 60 min — liczba rekordów z RT/ładunkiem zgodna z ETAPEM 02",
    nrow(kwalifikowane_60_panel) == oczekiwane_gradient("60 min", 3L, "Liczba_peptydow"),
    nrow(kwalifikowane_60_panel),
    oczekiwane_gradient("60 min", 3L, "Liczba_peptydow")
  ),
  kontrola(
    "Panel 60 min — jakość PEP + fragmenty zgodna z ETAPEM 02",
    nrow(jakosc_60_panel) == oczekiwane_gradient("60 min", 5L, "Liczba_peptydow"),
    nrow(jakosc_60_panel),
    oczekiwane_gradient("60 min", 5L, "Liczba_peptydow")
  ),
  kontrola(
    "Top 2 rozszerzona jest podzbiorem Top 3 rozszerzonej",
    nrow(top2_nie_w_top3) == 0,
    nrow(top2_nie_w_top3),
    0
  ),
  kontrola(
    "Top 2 — żaden przewidywany peptyd nie był obserwowany w 30 min",
    nrow(predykcja_obecna_w_30_top2) == 0,
    nrow(predykcja_obecna_w_30_top2),
    0
  ),
  kontrola(
    "Top 3 — żaden przewidywany peptyd nie był obserwowany w 30 min",
    nrow(predykcja_obecna_w_30_top3) == 0,
    nrow(predykcja_obecna_w_30_top3),
    0
  ),
  kontrola(
    "Wszystkie rekordy jakości 30 min mają skończony średni Fr.Score top3",
    all(is.finite(jakosc_30_panel$Max_sredni_score_top3)),
    sum(!is.finite(jakosc_30_panel$Max_sredni_score_top3)),
    0
  ),
  kontrola(
    "Wszystkie rekordy jakości 60 min mają skończony średni Fr.Score top3",
    all(is.finite(jakosc_60_panel$Max_sredni_score_top3)),
    sum(!is.finite(jakosc_60_panel$Max_sredni_score_top3)),
    0
  ),
  sprawdz_wariant(wariant_top2),
  sprawdz_wariant(wariant_top3)
)

snapshot_wynikow <- tibble::tribble(
  ~Miara, ~Wartosc,
  "Panel_obserwowane_30_przed_jakoscia", nrow(obserwowane_30_panel),
  "Panel_jakosc_30", nrow(jakosc_30_panel),
  "Panel_jakosc_60", nrow(jakosc_60_panel),
  "Panel_prawdziwie_tylko_60", nrow(
    jakosc_60_panel %>%
      anti_join(obserwowane_30_panel %>% select(Bialko, Peptyd), by = c("Bialko", "Peptyd"))
  ),
  "Walidacja_iRT_FASTA_N", nrow(wspolne_RT),
  "Walidacja_iRT_PANEL_N", sum(wspolne_RT$W_panelu_741),
  "Top2_dodane_RT30_obserwowany", nrow(wariant_top2$dodane_obserwowany),
  "Top2_dodane_RT30_przewidywany", nrow(wariant_top2$dodane_predykcja),
  "Top2_rozszerzona_peptydy", nrow(wariant_top2$rozszerzona),
  "Top2_rozszerzona_bialka", n_distinct(wariant_top2$rozszerzona$Bialko),
  "Top3_dodane_RT30_obserwowany", nrow(wariant_top3$dodane_obserwowany),
  "Top3_dodane_RT30_przewidywany", nrow(wariant_top3$dodane_predykcja),
  "Top3_rozszerzona_peptydy", nrow(wariant_top3$rozszerzona),
  "Top3_rozszerzona_bialka", n_distinct(wariant_top3$rozszerzona$Bialko)
)

zapisz_excel_wiele(
  list(
    Kontrole = kontrole,
    Snapshot_wynikow = snapshot_wynikow,
    Top2_nie_w_Top3 = top2_nie_w_top3,
    Predykcja_obecna_30_Top2 = predykcja_obecna_w_30_top2,
    Predykcja_obecna_30_Top3 = predykcja_obecna_w_30_top3
  ),
  file.path(folder_kontrola, "KONTROLA_ETAP_05_FINAL.xlsx")
)

readr::write_excel_csv2(
  kontrole,
  file.path(folder_kontrola, "KONTROLA_ETAP_05_FINAL.csv")
)

if (any(kontrole$Krytyczna & !kontrole$Wynik)) {
  stop(
    paste(
      "Co najmniej jedna krytyczna kontrola ETAPU 05 zakończyła się FALSE.",
      "Sprawdź folder 04_Kontrola."
    )
  )
}


# ------------------------------------------------------------------------------
# 11. EKSPORT TABEL
# ------------------------------------------------------------------------------

zapisz_excel_wiele(
  list(
    Parametry_kalibracji = kalibracja_opis,
    Wzory = wzory,
    Podsumowanie_bledu = podsumowanie_bledu,
    Wspolne_PANEL = wspolne_RT %>% filter(W_panelu_741)
  ),
  file.path(
    folder_tabele,
    "Tabela_01_Kalibracja_iRT_i_walidacja_RT30.xlsx"
  )
)

zapisz_excel_wiele(
  list(
    Top2_bazowa = wariant_top2$baza,
    Top2_rozszerzona = wariant_top2$rozszerzona,
    Top3_bazowa = wariant_top3$baza,
    Top3_rozszerzona = wariant_top3$rozszerzona
  ),
  file.path(
    folder_tabele,
    "Tabela_02_Listy_bazowe_i_rozszerzone_Top2_Top3.xlsx"
  )
)

# Alias zgodnościowy dla starszych skryptów/rycin.
file.copy(
  file.path(folder_tabele, "Tabela_02_Listy_bazowe_i_rozszerzone_Top2_Top3.xlsx"),
  file.path(folder_tabele, "Tabela_02_Listy_konserwatywne_i_rozszerzone_Top2_Top3.xlsx"),
  overwrite = TRUE
)

zapisz_excel_wiele(
  list(
    Top2_RT30_obs_jakosc60 = wariant_top2$dodane_obserwowany,
    Top2_prawdziwie_tylko60 = wariant_top2$dodane_predykcja,
    Top3_RT30_obs_jakosc60 = wariant_top3$dodane_obserwowany,
    Top3_prawdziwie_tylko60 = wariant_top3$dodane_predykcja
  ),
  file.path(
    folder_tabele,
    "Tabela_03_Peptydy_dodane_w_ETAPIE_05.xlsx"
  )
)

zapisz_excel_wiele(
  list(
    Top2_pula_RT30_obs = wariant_top2$pula_obserwowany,
    Top2_pula_tylko60 = wariant_top2$pula_tylko60,
    Top3_pula_RT30_obs = wariant_top3$pula_obserwowany,
    Top3_pula_tylko60 = wariant_top3$pula_tylko60
  ),
  file.path(
    folder_tabele,
    "Tabela_04_Pelne_pule_kandydatow_uzupelniajacych.xlsx"
  )
)

zapisz_excel_wiele(
  list(
    Status_Top2 = wariant_top2$status,
    Status_Top3 = wariant_top3$status,
    Podsumowanie_statusow = podsumowanie_statusow
  ),
  file.path(
    folder_tabele,
    "Tabela_05_Status_wszystkich_741_bialek_ETAP05.xlsx"
  )
)

zapisz_excel_wiele(
  list(
    Podsumowanie_list = podsumowanie_list,
    Dystrybucja = dystrybucja_peptydow,
    Podsumowanie_statusow = podsumowanie_statusow,
    Podsumowanie_bledu = podsumowanie_bledu
  ),
  file.path(
    folder_tabele,
    "Tabela_06_Podsumowanie_ETAP_05_FINAL.xlsx"
  )
)

zapisz_csv2(
  wariant_top2$baza,
  "ETAP_05_PANEL_Top2_lista_bazowa.csv"
)

zapisz_csv2(
  wariant_top2$rozszerzona,
  "ETAP_05_PANEL_Top2_lista_rozszerzona.csv"
)

zapisz_csv2(
  wariant_top2$dodane_obserwowany,
  "ETAP_05_PANEL_Top2_dodane_RT30_obserwowany_jakosc60.csv"
)

zapisz_csv2(
  wariant_top2$dodane_predykcja,
  "ETAP_05_PANEL_Top2_dodane_prawdziwie_tylko60_predykcja_iRT.csv"
)

zapisz_csv2(
  wariant_top3$baza,
  "ETAP_05_PANEL_Top3_lista_bazowa.csv"
)

# Aliasy zgodnościowe dla starszych skryptów rycin.
zapisz_csv2(
  wariant_top2$baza,
  "ETAP_05_PANEL_Top2_lista_konserwatywna.csv"
)
zapisz_csv2(
  wariant_top3$baza,
  "ETAP_05_PANEL_Top3_lista_konserwatywna.csv"
)

zapisz_csv2(
  wariant_top3$rozszerzona,
  "ETAP_05_PANEL_Top3_lista_rozszerzona.csv"
)

zapisz_csv2(
  wariant_top3$dodane_obserwowany,
  "ETAP_05_PANEL_Top3_dodane_RT30_obserwowany_jakosc60.csv"
)

zapisz_csv2(
  wariant_top3$dodane_predykcja,
  "ETAP_05_PANEL_Top3_dodane_prawdziwie_tylko60_predykcja_iRT.csv"
)

zapisz_csv2(
  bind_rows(wariant_top2$status, wariant_top3$status),
  "ETAP_05_status_741_bialek_Top2_Top3.csv"
)

zapisz_csv2(
  podsumowanie_list,
  "ETAP_05_podsumowanie_list.csv"
)

zapisz_csv2(
  podsumowanie_bledu,
  "ETAP_05_walidacja_iRT_podsumowanie.csv"
)


# ------------------------------------------------------------------------------
# 12. WYKRESY
# ------------------------------------------------------------------------------

p1 <- wspolne_RT %>%
  ggplot(
    aes(
      x = RT30_obserwowany,
      y = RT30_przewidziany,
      colour = W_panelu_741
    )
  ) +
  geom_point(alpha = 0.18, size = 0.8) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.8) +
  scale_colour_discrete(
    labels = c("Poza panelem", "Panel 741"),
    name = "Zakres"
  ) +
  labs(
    title = "Walidacja przewidywania RT30 z RT60 przez iRT",
    subtitle = paste(
      "Peptydy human-only spełniające PEP <= 0,01",
      "i co najmniej 3 kwalifikujące się fragmenty w obu gradientach"
    ),
    x = "RT30 obserwowany [min]",
    y = "RT30 przewidywany [min]",
    caption = "Linia y = x oznacza idealną zgodność."
  ) +
  theme_prm()

zapisz_wykres(
  p1,
  "Wykres_01_RT30_obserwowany_vs_przewidywany",
  9,
  6.5
)

p2 <- wspolne_RT %>%
  ggplot(aes(x = Blad_predykcji_min, fill = W_panelu_741)) +
  geom_histogram(bins = 60, alpha = 0.75, position = "identity") +
  geom_vline(xintercept = 0, linewidth = 0.8) +
  facet_wrap(~W_panelu_741, scales = "free_y") +
  labs(
    title = "Rozkład błędu predykcji RT30",
    subtitle = "Błąd = RT30 przewidywany - RT30 obserwowany",
    x = "Błąd [min]",
    y = "Liczba peptydów",
    fill = "Panel 741"
  ) +
  theme_prm()

zapisz_wykres(
  p2,
  "Wykres_02_Rozklad_bledu_predykcji_RT30",
  10,
  6
)

zrodla_wykres <- bind_rows(
  wariant_top2$rozszerzona,
  wariant_top3$rozszerzona
) %>%
  mutate(
    Kategoria = case_when(
      Typ_rekordu == "bazowy_ETAP03" ~
        "Bazowy: RT30 i jakość 30 min",
      Typ_rekordu == "dodany_RT30_obserwowany_jakosc60" ~
        "Dodany: RT30 obserwowany, jakość 60 min",
      TRUE ~
        "Dodany: prawdziwie tylko 60 min, RT30 przewidywany"
    )
  ) %>%
  count(Typ_listy, Kategoria, name = "Liczba_peptydow")

p3 <- zrodla_wykres %>%
  ggplot(
    aes(
      x = Typ_listy,
      y = Liczba_peptydow,
      fill = Kategoria
    )
  ) +
  geom_col(width = 0.65) +
  geom_text(
    aes(
      label = if_else(
        Liczba_peptydow > 0,
        scales::label_number(big.mark = " ")(Liczba_peptydow),
        ""
      )
    ),
    position = position_stack(vjust = 0.5),
    size = 3.5
  ) +
  scale_y_continuous(
    labels = scales::label_number(big.mark = " ")
  ) +
  labs(
    title = "Źródło RT i jakości w rozszerzonych listach panelowych",
    x = NULL,
    y = "Liczba peptydów",
    fill = "Typ rekordu"
  ) +
  theme_prm()

zapisz_wykres(
  p3,
  "Wykres_03_Zrodlo_RT_i_jakosci_Top2_Top3",
  11,
  6.5
)

p4 <- podsumowanie_list %>%
  select(
    Typ_listy,
    Wariant,
    Liczba_bialek,
    Liczba_peptydow
  ) %>%
  pivot_longer(
    cols = c(Liczba_bialek, Liczba_peptydow),
    names_to = "Miara",
    values_to = "Liczba"
  ) %>%
  mutate(
    Miara = recode(
      Miara,
      Liczba_bialek = "Białka",
      Liczba_peptydow = "Peptydy"
    )
  ) %>%
  ggplot(
    aes(
      x = Wariant,
      y = Liczba,
      fill = Miara
    )
  ) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65
  ) +
  geom_text(
    aes(label = scales::label_number(big.mark = " ")(Liczba)),
    position = position_dodge(width = 0.75),
    vjust = -0.25,
    size = 3.5
  ) +
  facet_wrap(~Typ_listy, scales = "free_y") +
  scale_y_continuous(
    labels = scales::label_number(big.mark = " "),
    expand = expansion(mult = c(0, 0.16))
  ) +
  labs(
    title = "Wpływ rozszerzenia na liczbę białek i peptydów",
    x = NULL,
    y = "Liczba",
    fill = "Miara"
  ) +
  theme_prm() +
  theme(axis.text.x = element_text(angle = 8, hjust = 0.5))

zapisz_wykres(
  p4,
  "Wykres_04_Bazowa_vs_rozszerzona",
  11,
  6.5
)

p5 <- dystrybucja_peptydow %>%
  mutate(
    Liczba_peptydow_na_bialko = factor(
      Liczba_peptydow_na_bialko,
      levels = 1:3
    ),
    Etykieta = paste0(Liczba_bialek, "\n", Procent_bialek, "%")
  ) %>%
  ggplot(
    aes(
      x = Liczba_peptydow_na_bialko,
      y = Liczba_bialek,
      fill = Typ_listy
    )
  ) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = Etykieta), vjust = -0.25, size = 4) +
  facet_wrap(~Typ_listy, scales = "free_y") +
  scale_y_continuous(
    labels = scales::label_number(big.mark = " "),
    expand = expansion(mult = c(0, 0.20))
  ) +
  labs(
    title = "Liczba peptydów na białko po rozszerzeniu",
    x = "Liczba peptydów na białko",
    y = "Liczba białek"
  ) +
  theme_prm()

zapisz_wykres(
  p5,
  "Wykres_05_Dystrybucja_peptydow_na_bialko",
  10,
  6
)


# ------------------------------------------------------------------------------
# 13. INDEKS I RAPORT
# ------------------------------------------------------------------------------

indeks_plikow <- tibble::tribble(
  ~Folder, ~Plik, ~Opis,
  "01_Tabele_Excel",
  "Tabela_01_Kalibracja_iRT_i_walidacja_RT30.xlsx",
  "Parametry kalibracji, wzory i walidacja przewidywania RT30.",
  "01_Tabele_Excel",
  "Tabela_02_Listy_bazowe_i_rozszerzone_Top2_Top3.xlsx",
  "Główne listy Top 2 i Top 3 przed i po rozszerzeniu.",
  "01_Tabele_Excel",
  "Tabela_03_Peptydy_dodane_w_ETAPIE_05.xlsx",
  "Rozdzielenie dodatków z obserwowanym RT30 i prawdziwych tylko-60.",
  "01_Tabele_Excel",
  "Tabela_04_Pelne_pule_kandydatow_uzupelniajacych.xlsx",
  "Pełne pule rankingowe kandydatów dla Top 2 i Top 3.",
  "01_Tabele_Excel",
  "Tabela_05_Status_wszystkich_741_bialek_ETAP05.xlsx",
  "Dokładnie jeden status dla każdego białka i każdego wariantu listy.",
  "01_Tabele_Excel",
  "Tabela_06_Podsumowanie_ETAP_05_FINAL.xlsx",
  "Podsumowania list, błędu predykcji, statusów i dystrybucji.",
  "02_CSV",
  "ETAP_05_PANEL_Top2_lista_rozszerzona.csv",
  "Rozszerzona lista panelowa Top 2.",
  "02_CSV",
  "ETAP_05_PANEL_Top3_lista_rozszerzona.csv",
  "Rozszerzona lista panelowa Top 3.",
  "04_Kontrola",
  "KONTROLA_ETAP_05_FINAL.xlsx",
  "Pełny audyt krytycznych reguł ETAPU 05."
)

zapisz_excel(
  indeks_plikow,
  file.path(folder_wyniki, "00_INDEKS_PLIKOW_ETAP_05_FINAL.xlsx")
)

raport <- c(
  "ETAP 05 — KALKULATOR iRT I LISTY ROZSZERZONE — EXPORT-QUANT — FINAL HUMAN-ONLY",
  paste("Data wykonania:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "NAJWAŻNIEJSZA KOREKTA:",
  paste(
    "Peptyd jest uznawany za prawdziwie tylko-60 min dopiero po anti-join",
    "względem wszystkich peptydów zaobserwowanych w 30 min,",
    "przed filtrami PEP i fragmentacyjnymi."
  ),
  "",
  "WALIDACJA iRT:",
  paste(
    "FASTA human-only:",
    nrow(wspolne_RT),
    "wspólnych peptydów;"
  ),
  paste(
    "panel 741:",
    sum(wspolne_RT$W_panelu_741),
    "wspólnych peptydów."
  ),
  paste(
    "P95 błędu użyty orientacyjnie:",
    round(margines_95_uzyty, 6),
    "min."
  ),
  "",
  "TOP 2:",
  paste(
    "Baza:",
    nrow(wariant_top2$baza),
    "peptydów /",
    n_distinct(wariant_top2$baza$Bialko),
    "białek."
  ),
  paste(
    "Dodane z obserwowanym RT30 i jakością 60 min:",
    nrow(wariant_top2$dodane_obserwowany)
  ),
  paste(
    "Dodane prawdziwie tylko-60 z predykcją RT30:",
    nrow(wariant_top2$dodane_predykcja)
  ),
  paste(
    "Rozszerzona:",
    nrow(wariant_top2$rozszerzona),
    "peptydów /",
    n_distinct(wariant_top2$rozszerzona$Bialko),
    "białek."
  ),
  "",
  "TOP 3:",
  paste(
    "Baza:",
    nrow(wariant_top3$baza),
    "peptydów /",
    n_distinct(wariant_top3$baza$Bialko),
    "białek."
  ),
  paste(
    "Dodane z obserwowanym RT30 i jakością 60 min:",
    nrow(wariant_top3$dodane_obserwowany)
  ),
  paste(
    "Dodane prawdziwie tylko-60 z predykcją RT30:",
    nrow(wariant_top3$dodane_predykcja)
  ),
  paste(
    "Rozszerzona:",
    nrow(wariant_top3$rozszerzona),
    "peptydów /",
    n_distinct(wariant_top3$rozszerzona$Bialko),
    "białek."
  ),
  "",
  paste(
    "Wszystkie krytyczne kontrole:",
    ifelse(all(kontrole$Wynik[kontrole$Krytyczna]), "TRUE", "FALSE")
  ),
  "",
  paste(
    "UWAGA: orientacyjny zakres +/- P95 nie jest finalnym oknem scheduled PRM.",
    "Listy wymagają kontroli w Skyline."
  )
)

writeLines(
  raport,
  file.path(folder_kontrola, "RAPORT_ETAP_05_FINAL.txt"),
  useBytes = TRUE
)

cat("\nETAP 05 zakończony poprawnie.\n")
cat("Folder wynikowy:\n", folder_wyniki, "\n\n")
cat(
  "Top 2 rozszerzona:",
  nrow(wariant_top2$rozszerzona),
  "peptydów /",
  n_distinct(wariant_top2$rozszerzona$Bialko),
  "białek.\n"
)
cat(
  "Top 3 rozszerzona:",
  nrow(wariant_top3$rozszerzona),
  "peptydów /",
  n_distinct(wariant_top3$rozszerzona$Bialko),
  "białek.\n"
)
cat(
  "Prawdziwie przewidywane RT30:",
  nrow(wariant_top2$dodane_predykcja),
  "(Top 2) i",
  nrow(wariant_top3$dodane_predykcja),
  "(Top 3).\n"
)

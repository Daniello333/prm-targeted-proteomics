# ==============================================================================
# 01. PRZYGOTOWANIE FASTA I TRAWIENIE IN SILICO
#
# Praca magisterska:
# „Planowanie celowanych eksperymentów proteomicznych z wykorzystaniem środowiska R”
# Autor: Daniel Drabik
#
# Skrypt przygotowuje część sekwencyjną analizy. Wczytuje referencyjny proteom
# Homo sapiens oraz panel 741 białek, wykonuje trawienie in silico trypsyną i
# Lys-C dla 0 i 1 DPMC, a następnie wyznacza peptydy unikatowe względem całej
# analizowanej FASTA. Dla wariantu Trypsyna 1 DPMC stosowany jest dodatkowo filtr
# potencjalnych peptydów ragged.
#
# Wyniki tego skryptu są wejściem do etapu 02. Katalog wynikowy jest tworzony
# automatycznie i otrzymuje znacznik daty oraz czasu uruchomienia.
# ==============================================================================

# ----------------------------------------------------------------------------
# 0. Pakiety
# ----------------------------------------------------------------------------

pakiety <- c("tidyverse", "readxl", "openxlsx")

for (p in pakiety) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}

library(tidyverse)
library(readxl)
library(openxlsx)

# ----------------------------------------------------------------------------
# 1. Parametry analizy i katalogi wynikowe
# ----------------------------------------------------------------------------

# Można podać pełne ścieżki ręcznie. NULL uruchamia bezpieczne wyszukiwanie
# najbardziej aktualnego zgodnego pliku w bieżącym katalogu i podkatalogach.
plik_fasta <- NULL
plik_panel <- NULL

oczekiwana_liczba_rekordow_Homo_sapiens <- 20663L
oczekiwana_liczba_standardow_iRT <- 1L
oczekiwana_liczba_bialek_panelu <- 741L

data_fasta_z_nazwy_pliku <- "2025-11-19"
zrodlo_fasta <- paste(
  "UniProt reference proteome UP000005640, Homo sapiens (taxon 9606);",
  "plik zawiera dodatkowo techniczny rekord Biognosys iRT"
)

min_dlugosc <- 7L
max_dlugosc <- 25L
max_dlugosc_kontrolna <- 35L

# Filtr ragged jest stosowany wyłącznie do wariantu
# Trypsyna 1 DPMC.
usun_ragged_trypsyna_1DPMC <- TRUE

# Sekwencje uznawane za potencjalne ragged ends:
# KK, KR, RK lub RR na N- lub C-końcu peptydu.
wzorzec_ragged_N <- "^(KK|KR|RK|RR)"
wzorzec_ragged_C <- "(KK|KR|RK|RR)$"

znacznik <- format(Sys.time(), "%Y%m%d_%H%M%S")
folder_wynikowy <- paste0(
  "WYNIKI_ETAP_00_01_KOREKTY_SKK_HUMAN_ONLY_",
  znacznik
)
dir.create(folder_wynikowy, recursive = TRUE, showWarnings = FALSE)

folder_csv <- file.path(folder_wynikowy, "01_CSV")
folder_xlsx <- file.path(folder_wynikowy, "02_XLSX")
folder_kontrola <- file.path(folder_wynikowy, "03_KONTROLA")

dir.create(folder_csv, recursive = TRUE, showWarnings = FALSE)
dir.create(folder_xlsx, recursive = TRUE, showWarnings = FALSE)
dir.create(folder_kontrola, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------------
# 2. Funkcje pomocnicze
# ----------------------------------------------------------------------------

# Sprawdza dostępność pliku wejściowego i przerywa analizę, jeżeli plik nie istnieje.
sprawdz_plik <- function(sciezka) {
  if (!file.exists(sciezka)) {
    stop("Nie znaleziono pliku: ", sciezka)
  }
}

# Wybiera najnowszy istniejący plik na podstawie czasu modyfikacji.
wybierz_najnowszy_plik <- function(sciezki) {
  sciezki <- unique(sciezki[file.exists(sciezki)])
  if (length(sciezki) == 0) return(NULL)
  sciezki[which.max(file.info(sciezki)$mtime)]
}

# Wybiera plik wskazany ręcznie albo automatycznie odnajduje najnowszy plik zgodny z podanym wzorcem.
znajdz_plik_wejsciowy <- function(sciezka_reczna, wzorce, opis) {
  if (!is.null(sciezka_reczna)) {
    sprawdz_plik(sciezka_reczna)
    return(normalizePath(sciezka_reczna, winslash = "/", mustWork = TRUE))
  }

  wszystkie <- list.files(
    path = ".",
    recursive = TRUE,
    full.names = TRUE,
    include.dirs = FALSE
  )

  dopasowane <- wszystkie[
    vapply(
      basename(wszystkie),
      function(nazwa) any(stringr::str_detect(nazwa, wzorce)),
      logical(1)
    )
  ]

  wybrany <- wybierz_najnowszy_plik(dopasowane)
  if (is.null(wybrany)) {
    stop("Nie znaleziono pliku wejściowego: ", opis)
  }

  normalizePath(wybrany, winslash = "/", mustWork = TRUE)
}

# Zapisuje tabelę jako plik CSV2 czytelny w polskiej konfiguracji programu Excel.
zapisz_csv2 <- function(df, nazwa_pliku, folder_docelowy = folder_csv) {
  readr::write_excel_csv2(
    df,
    file.path(folder_docelowy, nazwa_pliku),
    na = ""
  )
}

# Sprawdza, czy sekwencja zawiera wyłącznie standardowe jednoliterowe oznaczenia aminokwasów.
czy_standardowa_sekwencja <- function(peptyd) {
  stringr::str_detect(peptyd, "^[ACDEFGHIKLMNPQRSTVWY]+$")
}

# Sprawdza obecność motywu ragged na końcu N peptydu.
czy_ragged_N <- function(peptyd) {
  stringr::str_detect(peptyd, wzorzec_ragged_N)
}

# Sprawdza obecność motywu ragged na końcu C peptydu.
czy_ragged_C <- function(peptyd) {
  stringr::str_detect(peptyd, wzorzec_ragged_C)
}

# Oznacza sekwencje z motywami KK, KR, RK lub RR na jednym z końców peptydu.
czy_ragged <- function(peptyd) {
  czy_ragged_N(peptyd) | czy_ragged_C(peptyd)
}

# Przypisuje peptyd do kategorii długości wykorzystywanej w tabelach podsumowujących.
kategoria_dlugosci <- function(dlugosc) {
  case_when(
    between(dlugosc, 7L, 10L) ~ "7-10 aa",
    between(dlugosc, 11L, 15L) ~ "11-15 aa",
    between(dlugosc, 16L, 20L) ~ "16-20 aa",
    between(dlugosc, 21L, 25L) ~ "21-25 aa",
    between(dlugosc, 26L, 30L) ~ "26-30 aa",
    between(dlugosc, 31L, 35L) ~ "31-35 aa",
    TRUE ~ "poza analizowanym zakresem"
  )
}

# ----------------------------------------------------------------------------
# 3. Obliczanie mas monoizotopowych i teoretycznych wartości m/z
# ----------------------------------------------------------------------------

masa_reszt <- c(
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

masa_H2O <- 18.010564684
masa_protonu <- 1.007276466812

# Oblicza neutralną masę monoizotopową peptydu na podstawie sumy mas reszt aminokwasowych i cząsteczki wody.
oblicz_mase_peptydu <- function(peptyd) {
  aa <- strsplit(peptyd, "", fixed = TRUE)[[1]]
  if (any(!aa %in% names(masa_reszt))) return(NA_real_)
  sum(masa_reszt[aa]) + masa_H2O
}

# Przelicza neutralną masę peptydu na wartość m/z dla wskazanego stanu naładowania.
oblicz_mz <- function(masa_neutralna, z) {
  (masa_neutralna + z * masa_protonu) / z
}

# Dodaje do tabeli neutralną masę monoizotopową oraz teoretyczne wartości m/z dla kilku stanów naładowania.
dodaj_mz_teoretyczne <- function(df) {
  df %>%
    mutate(
      Masa_monoizotopowa = map_dbl(Peptyd, oblicz_mase_peptydu),
      `m/z z=2` = round(oblicz_mz(Masa_monoizotopowa, 2), 6),
      `m/z z=3` = round(oblicz_mz(Masa_monoizotopowa, 3), 6),
      `m/z z=4` = round(oblicz_mz(Masa_monoizotopowa, 4), 6),
      `m/z z=5` = round(oblicz_mz(Masa_monoizotopowa, 5), 6)
    )
}

# ----------------------------------------------------------------------------
# 4. Wczytanie i kontrola pliku FASTA
# ----------------------------------------------------------------------------

# Wczytuje rekordy FASTA, porządkuje nagłówki i sekwencje oraz przygotowuje tabelę białek.
wczytaj_fasta <- function(sciezka) {
  sprawdz_plik(sciezka)

  linie <- readLines(sciezka, warn = FALSE)
  naglowki_idx <- grep("^>", linie)

  if (length(naglowki_idx) == 0) {
    stop("Plik FASTA nie zawiera naglowkow zaczynajacych sie od '>'.")
  }

  map_dfr(seq_along(naglowki_idx), function(i) {
    start <- naglowki_idx[i]
    koniec <- if (i < length(naglowki_idx)) {
      naglowki_idx[i + 1] - 1L
    } else {
      length(linie)
    }

    naglowek <- linie[start]

    sekwencja <- if (start < koniec) {
      paste0(linie[(start + 1):koniec], collapse = "")
    } else {
      ""
    }

    accession <- str_match(naglowek, "^>(?:sp|tr)\\|([^|]+)\\|")[, 2]
    if (is.na(accession)) {
      accession <- str_match(naglowek, "^>(\\S+)")[, 2]
    }

    nazwa <- naglowek %>%
      str_remove("^>") %>%
      str_remove("^(sp|tr)\\|[^|]+\\|[^ ]+\\s*") %>%
      str_remove("\\sOS=.*$") %>%
      str_trim()

    if (is.na(nazwa) || nazwa == "") nazwa <- accession

    tibble(
      Bialko = accession,
      Nazwa_bialka = nazwa,
      Naglowek_FASTA = naglowek,
      Sekwencja = toupper(str_replace_all(sekwencja, "\\s+", ""))
    )
  }) %>%
    filter(!is.na(Bialko), Bialko != "", Sekwencja != "") %>%
    distinct(Bialko, .keep_all = TRUE)
}

# ----------------------------------------------------------------------------
# 5. Trawienie in silico
# ----------------------------------------------------------------------------

# Generuje produkty trawienia pojedynczej sekwencji zgodnie z regułą wybranego enzymu i liczbą DPMC.
trawienie_bialka <- function(sekwencja, enzym = "trypsyna", max_dpmc = 0L) {
  aa <- strsplit(sekwencja, "", fixed = TRUE)[[1]]
  n <- length(aa)

  if (n < 1) return(character(0))

  reszty_ciecia <- if (enzym == "trypsyna") c("K", "R") else "K"

  # Definicja odpowiada ustawieniu Trypsin: brak ciecia przed prolina.
  miejsca_ciecia <- integer(0)
  if (n >= 2) {
    miejsca_ciecia <- which(
      aa[-n] %in% reszty_ciecia & aa[-1] != "P"
    )
  }

  granice <- c(0L, miejsca_ciecia, n)
  liczba_fragmentow <- length(granice) - 1L
  peptydy <- character(0)

  for (dpmc in 0:max_dpmc) {
    if ((liczba_fragmentow - dpmc) <= 0) next

    for (i in seq_len(liczba_fragmentow - dpmc)) {
      poczatek <- granice[i] + 1L
      koniec <- granice[i + 1L + dpmc]
      peptydy <- c(peptydy, paste0(aa[poczatek:koniec], collapse = ""))
    }
  }

  unique(peptydy)
}

# Wykonuje trawienie wszystkich białek i zapisuje relacje pomiędzy peptydami a białkami źródłowymi.
generuj_przypisania_peptydow <- function(
    fasta_df,
    enzym,
    max_dpmc,
    min_len = min_dlugosc,
    max_len = max_dlugosc,
    usun_ragged = FALSE
) {
  fasta_df %>%
    mutate(
      Peptydy = map(
        Sekwencja,
        ~ trawienie_bialka(.x, enzym = enzym, max_dpmc = max_dpmc)
      )
    ) %>%
    select(Bialko, Nazwa_bialka, Peptydy) %>%
    unnest(Peptydy) %>%
    rename(Peptyd = Peptydy) %>%
    distinct(Bialko, Nazwa_bialka, Peptyd) %>%
    mutate(
      Dlugosc_peptydu = nchar(Peptyd),
      Standardowe_aminokwasy = czy_standardowa_sekwencja(Peptyd),
      Zawiera_metionine = str_detect(Peptyd, "M"),
      Ragged_N = czy_ragged_N(Peptyd),
      Ragged_C = czy_ragged_C(Peptyd),
      Ragged = Ragged_N | Ragged_C,
      Kategoria_dlugosci = kategoria_dlugosci(Dlugosc_peptydu)
    ) %>%
    filter(
      Dlugosc_peptydu >= min_len,
      Dlugosc_peptydu <= max_len,
      Standardowe_aminokwasy,
      !Zawiera_metionine
    ) %>%
    {
      if (usun_ragged) filter(., !Ragged) else .
    }
}

# Pozostawia peptydy przypisane do dokładnie jednego białka w całym proteomie.
wybierz_unikatowe <- function(przypisania_df) {
  przypisania_df %>%
    group_by(Peptyd) %>%
    mutate(Liczba_bialek_dla_sekwencji = n_distinct(Bialko)) %>%
    ungroup() %>%
    filter(Liczba_bialek_dla_sekwencji == 1L) %>%
    select(-Liczba_bialek_dla_sekwencji)
}

# ----------------------------------------------------------------------------
# 6. Funkcje podsumowujące
# ----------------------------------------------------------------------------

# Oblicza podstawowe liczebności dla jednego wariantu trawienia.
podsumuj_wariant <- function(
    unikatowe_df,
    fasta_df,
    nazwa_wariantu,
    panel_ids
) {
  liczba_bialek_fasta <- n_distinct(fasta_df$Bialko)
  liczba_bialek_z_peptydem <- n_distinct(unikatowe_df$Bialko)

  panel_df <- unikatowe_df %>% filter(Bialko %in% panel_ids)

  tibble(
    Wariant = nazwa_wariantu,
    Bialka_w_FASTA = liczba_bialek_fasta,
    Bialka_z_min_1_peptydem_unikatowym = liczba_bialek_z_peptydem,
    Bialka_bez_peptydu_unikatowego = liczba_bialek_fasta - liczba_bialek_z_peptydem,
    Peptydy_unikatowe = nrow(unikatowe_df),
    Pokrycie_FASTA_proc = round(100 * liczba_bialek_z_peptydem / liczba_bialek_fasta, 4),
    Panel_bialka_z_min_1_peptydem = n_distinct(panel_df$Bialko),
    Panel_bialka_bez_peptydu = length(setdiff(panel_ids, unique(panel_df$Bialko))),
    Panel_peptydy_unikatowe = nrow(panel_df),
    Pokrycie_panelu_proc = round(100 * n_distinct(panel_df$Bialko) / length(panel_ids), 4)
  )
}

# Zestawia liczbę peptydów w przyjętych kategoriach długości.
podsumuj_dlugosci <- function(unikatowe_df, nazwa_wariantu, zakres_analizy) {
  poziomy <- c("7-10 aa", "11-15 aa", "16-20 aa", "21-25 aa")

  unikatowe_df %>%
    mutate(
      Kategoria_dlugosci = factor(
        kategoria_dlugosci(nchar(Peptyd)),
        levels = poziomy
      )
    ) %>%
    filter(!is.na(Kategoria_dlugosci)) %>%
    count(Kategoria_dlugosci, name = "Liczba_peptydow") %>%
    complete(
      Kategoria_dlugosci = factor(poziomy, levels = poziomy),
      fill = list(Liczba_peptydow = 0L)
    ) %>%
    mutate(
      Wariant = nazwa_wariantu,
      Zakres_analizy = zakres_analizy,
      Udzial_proc = round(100 * Liczba_peptydow / sum(Liczba_peptydow), 4)
    ) %>%
    select(
      Zakres_analizy,
      Wariant,
      Kategoria_dlugosci,
      Liczba_peptydow,
      Udzial_proc
    )
}

# Wyznacza rozkład liczby peptydów unikatowych przypadających na białko.
podsumuj_peptydy_na_bialko <- function(unikatowe_df, nazwa_wariantu, zakres_analizy, wszystkie_bialka) {
  liczebnosci <- unikatowe_df %>%
    count(Bialko, name = "Liczba_peptydow_unikatowych")

  tibble(Bialko = wszystkie_bialka) %>%
    left_join(liczebnosci, by = "Bialko") %>%
    mutate(
      Liczba_peptydow_unikatowych = replace_na(Liczba_peptydow_unikatowych, 0L),
      Kategoria = case_when(
        Liczba_peptydow_unikatowych == 0L ~ "0",
        Liczba_peptydow_unikatowych == 1L ~ "1",
        Liczba_peptydow_unikatowych == 2L ~ "2",
        between(Liczba_peptydow_unikatowych, 3L, 5L) ~ "3-5",
        between(Liczba_peptydow_unikatowych, 6L, 10L) ~ "6-10",
        between(Liczba_peptydow_unikatowych, 11L, 20L) ~ "11-20",
        Liczba_peptydow_unikatowych > 20L ~ ">20"
      ),
      Kategoria = factor(Kategoria, levels = c("0", "1", "2", "3-5", "6-10", "11-20", ">20"))
    ) %>%
    count(Kategoria, name = "Liczba_bialek") %>%
    complete(
      Kategoria = factor(c("0", "1", "2", "3-5", "6-10", "11-20", ">20"),
                         levels = c("0", "1", "2", "3-5", "6-10", "11-20", ">20")),
      fill = list(Liczba_bialek = 0L)
    ) %>%
    mutate(
      Wariant = nazwa_wariantu,
      Zakres_analizy = zakres_analizy,
      Udzial_proc = round(100 * Liczba_bialek / sum(Liczba_bialek), 4)
    ) %>%
    select(Zakres_analizy, Wariant, Kategoria, Liczba_bialek, Udzial_proc)
}

# ----------------------------------------------------------------------------
# 7. Wczytanie danych wejściowych
# ----------------------------------------------------------------------------

plik_fasta <- znajdz_plik_wejsciowy(
  sciezka_reczna = plik_fasta,
  wzorce = c(
    "^2025_11_19_homo_sapiens_20663_UP000005640_9606_iRT(?:\\(1\\))?[.]fasta$"
  ),
  opis = "FASTA Homo sapiens UP000005640 z dołączonym standardem iRT"
)

plik_panel <- znajdz_plik_wejsciowy(
  sciezka_reczna = plik_panel,
  wzorce = c(
    "^ALL_candidates_741[.]xlsx$",
    "^ALL_candidates_741_REKONSTRUKCJA_ID[.]xlsx$"
  ),
  opis = "panel 741 białek"
)

fasta_wszystkie <- wczytaj_fasta(plik_fasta) %>%
  mutate(
    Czy_Homo_sapiens = str_detect(
      Naglowek_FASTA,
      "(^|\\s)OX=9606(\\s|$)"
    ),
    Czy_standard_iRT =
      Bialko == "Biognosys|iRT-Kit_WR_fusion" |
      str_detect(Naglowek_FASTA, regex("iRT[-_ ]?Kit", ignore_case = TRUE))
  )

fasta_df <- fasta_wszystkie %>%
  filter(Czy_Homo_sapiens, !Czy_standard_iRT) %>%
  select(-Czy_Homo_sapiens, -Czy_standard_iRT)

fasta_irt <- fasta_wszystkie %>%
  filter(Czy_standard_iRT) %>%
  select(-Czy_Homo_sapiens, -Czy_standard_iRT)

fasta_niesklasyfikowane <- fasta_wszystkie %>%
  filter(!Czy_Homo_sapiens, !Czy_standard_iRT)

if (nrow(fasta_niesklasyfikowane) > 0) {
  stop(
    "FASTA zawiera rekordy, których nie sklasyfikowano ani jako Homo sapiens, ",
    "ani jako standard iRT: ",
    paste(fasta_niesklasyfikowane$Bialko, collapse = ", ")
  )
}

if (n_distinct(fasta_df$Bialko) != oczekiwana_liczba_rekordow_Homo_sapiens) {
  stop(
    "Nieoczekiwana liczba rekordów Homo sapiens. Oczekiwano ",
    oczekiwana_liczba_rekordow_Homo_sapiens,
    ", otrzymano ",
    n_distinct(fasta_df$Bialko),
    "."
  )
}

if (n_distinct(fasta_irt$Bialko) != oczekiwana_liczba_standardow_iRT) {
  stop(
    "Nieoczekiwana liczba technicznych rekordów iRT. Oczekiwano ",
    oczekiwana_liczba_standardow_iRT,
    ", otrzymano ",
    n_distinct(fasta_irt$Bialko),
    "."
  )
}

# Kontrola dotyczy wyłącznie dokładnego rekordu technicznego.
# Nie wolno wyszukiwać ogólnego napisu "iRT" bez rozróżniania wielkości liter,
# ponieważ występuje on również w prawidłowych nazwach białek, np. SIRT1–SIRT7 i PIRT.
czy_techniczny_irt_w_human <-
  fasta_df$Bialko == "Biognosys|iRT-Kit_WR_fusion" |
  str_detect(
    fasta_df$Naglowek_FASTA,
    regex("^>Biognosys\\|iRT-Kit_WR_fusion(?:\\s|$)", ignore_case = FALSE)
  )

if (any(czy_techniczny_irt_w_human)) {
  stop("Techniczny rekord Biognosys iRT nadal znajduje się w zbiorze Homo sapiens.")
}

kontrola_rekordow_fasta <- tibble(
  Kategoria = c(
    "Wszystkie rekordy w pliku",
    "Rekordy Homo sapiens użyte w analizie",
    "Techniczne rekordy iRT wyłączone z analizy",
    "Rekordy niesklasyfikowane"
  ),
  Liczba_rekordow = c(
    n_distinct(fasta_wszystkie$Bialko),
    n_distinct(fasta_df$Bialko),
    n_distinct(fasta_irt$Bialko),
    n_distinct(fasta_niesklasyfikowane$Bialko)
  ),
  Uzyte_w_analizie_proteomu = c(FALSE, TRUE, FALSE, FALSE)
)

panel_surowy <- read_excel(plik_panel)
kolumna_id <- if ("ID" %in% colnames(panel_surowy)) {
  "ID"
} else if ("UniProt" %in% colnames(panel_surowy)) {
  "UniProt"
} else if ("Bialko" %in% colnames(panel_surowy)) {
  "Bialko"
} else {
  stop("W pliku panelowym nie znaleziono kolumny ID, UniProt ani Bialko.")
}

panel_ids <- panel_surowy %>%
  transmute(Bialko = str_trim(as.character(.data[[kolumna_id]]))) %>%
  filter(!is.na(Bialko), Bialko != "") %>%
  distinct() %>%
  pull(Bialko)

if (length(panel_ids) != oczekiwana_liczba_bialek_panelu) {
  stop(
    "Nieoczekiwana liczba identyfikatorów panelu. Oczekiwano ",
    oczekiwana_liczba_bialek_panelu,
    ", otrzymano ",
    length(panel_ids),
    "."
  )
}

panel_nieobecny_w_fasta <- setdiff(panel_ids, fasta_df$Bialko)

if (length(panel_nieobecny_w_fasta) > 0) {
  warning(
    "Nie wszystkie identyfikatory panelowe sa obecne w FASTA: ",
    paste(panel_nieobecny_w_fasta, collapse = ", ")
  )
}

cat("Plik FASTA:", plik_fasta, "\n")
cat("Wszystkie rekordy w pliku FASTA:", n_distinct(fasta_wszystkie$Bialko), "\n")
cat("Rekordy Homo sapiens użyte w analizie:", n_distinct(fasta_df$Bialko), "\n")
cat("Techniczne rekordy iRT wyłączone z analizy:", n_distinct(fasta_irt$Bialko), "\n")
cat("Liczba unikatowych identyfikatorow panelu:", length(panel_ids), "\n")
cat("Identyfikatory panelowe obecne w FASTA:", length(intersect(panel_ids, fasta_df$Bialko)), "\n")

# ----------------------------------------------------------------------------
# 8. Generowanie czterech wariantów trawienia i zastosowanie filtra ragged
# ----------------------------------------------------------------------------

warianty_def <- tribble(
  ~Klucz,                  ~Etykieta,                    ~Enzym,      ~DPMC, ~Usun_ragged,
  "Trypsyna_0MC",         "Trypsyna, 0 DPMC",           "trypsyna", 0L,    FALSE,
  "Trypsyna_1MC_PRZED",   "Trypsyna, 1 DPMC przed filtrem ragged", "trypsyna", 1L, FALSE,
  "Trypsyna_1MC",         "Trypsyna, 1 DPMC po filtrze ragged",    "trypsyna", 1L, usun_ragged_trypsyna_1DPMC,
  "LysC_0MC",             "LysC, 0 DPMC",               "lysc",     0L,    FALSE,
  "LysC_1MC",             "LysC, 1 DPMC",               "lysc",     1L,    FALSE
)

unikatowe_warianty <- list()
przypisania_warianty <- list()
podsumowania <- list()
podsumowania_dlugosci <- list()
podsumowania_dlugosci_panel <- list()
dystrybucje_bialek <- list()
dystrybucje_bialek_panel <- list()

for (i in seq_len(nrow(warianty_def))) {
  klucz <- warianty_def$Klucz[i]
  etykieta <- warianty_def$Etykieta[i]

  cat("\nGenerowanie wariantu:", etykieta, "\n")

  przypisania <- generuj_przypisania_peptydow(
    fasta_df = fasta_df,
    enzym = warianty_def$Enzym[i],
    max_dpmc = warianty_def$DPMC[i],
    min_len = min_dlugosc,
    max_len = max_dlugosc,
    usun_ragged = warianty_def$Usun_ragged[i]
  )

  unikatowe <- wybierz_unikatowe(przypisania)

  przypisania_warianty[[klucz]] <- przypisania
  unikatowe_warianty[[klucz]] <- unikatowe

  podsumowania[[klucz]] <- podsumuj_wariant(
    unikatowe,
    fasta_df,
    etykieta,
    panel_ids
  )

  podsumowania_dlugosci[[klucz]] <- podsumuj_dlugosci(
    unikatowe,
    etykieta,
    "Cała FASTA"
  )

  podsumowania_dlugosci_panel[[klucz]] <- podsumuj_dlugosci(
    unikatowe %>% filter(Bialko %in% panel_ids),
    etykieta,
    "Panel 741 białek"
  )

  dystrybucje_bialek[[klucz]] <- podsumuj_peptydy_na_bialko(
    unikatowe,
    etykieta,
    "Cała FASTA",
    fasta_df$Bialko
  )

  dystrybucje_bialek_panel[[klucz]] <- podsumuj_peptydy_na_bialko(
    unikatowe %>% filter(Bialko %in% panel_ids),
    etykieta,
    "Panel 741 białek",
    panel_ids
  )

  unikatowe_z_mz <- unikatowe %>%
    dodaj_mz_teoretyczne() %>%
    select(
      Bialko,
      Nazwa_bialka,
      Peptyd,
      Dlugosc_peptydu,
      Kategoria_dlugosci,
      Ragged_N,
      Ragged_C,
      Ragged,
      `m/z z=2`,
      `m/z z=3`,
      `m/z z=4`,
      `m/z z=5`
    ) %>%
    arrange(Bialko, Peptyd)

  # Nazwy zgodne z dalszymi etapami. W przypadku Trypsyny 1 MC plik zawiera
  # juz wynik po filtrze ragged.
  if (klucz != "Trypsyna_1MC_PRZED") {
    zapisz_csv2(
      unikatowe_z_mz,
      paste0("E1_peptydy_unikatowe_", klucz, ".csv"),
      folder_docelowy = folder_wynikowy
    )
  } else {
    zapisz_csv2(
      unikatowe_z_mz,
      "E1_peptydy_unikatowe_Trypsyna_1MC_PRZED_RAGGED.csv",
      folder_docelowy = folder_wynikowy
    )
  }
}

# Tabela kontrolna zawiera pięć wierszy, ponieważ pokazuje Trypsynę 1 DPMC
# zarówno przed, jak i po filtrze ragged.
podsumowanie_wariantow_z_kontrola <- bind_rows(podsumowania)

# Plik zgodnościowy „E1_podsumowanie_4_warianty.csv” musi zawierać dokładnie
# cztery końcowe warianty używane w dalszej analizie.
podsumowanie_wariantow <- bind_rows(
  podsumowania[c("Trypsyna_0MC", "Trypsyna_1MC", "LysC_0MC", "LysC_1MC")]
)

podsumowanie_dlugosci_fasta <- bind_rows(podsumowania_dlugosci)
podsumowanie_dlugosci_panel <- bind_rows(podsumowania_dlugosci_panel)
dystrybucja_bialek_fasta <- bind_rows(dystrybucje_bialek)
dystrybucja_bialek_panel <- bind_rows(dystrybucje_bialek_panel)

# ----------------------------------------------------------------------------
# 9. Porównanie wariantu Trypsyna 1 DPMC przed i po filtrze ragged
# ----------------------------------------------------------------------------

trypsyna_1_przed <- unikatowe_warianty[["Trypsyna_1MC_PRZED"]]
trypsyna_1_po <- unikatowe_warianty[["Trypsyna_1MC"]]

usuniete_ragged <- trypsyna_1_przed %>%
  anti_join(
    trypsyna_1_po %>% select(Bialko, Peptyd),
    by = c("Bialko", "Peptyd")
  ) %>%
  mutate(
    Dlugosc_peptydu = nchar(Peptyd),
    Ragged_N = czy_ragged_N(Peptyd),
    Ragged_C = czy_ragged_C(Peptyd),
    Ragged = Ragged_N | Ragged_C,
    Motyw_N = if_else(Ragged_N, str_sub(Peptyd, 1, 2), NA_character_),
    Motyw_C = if_else(Ragged_C, str_sub(Peptyd, -2, -1), NA_character_)
  ) %>%
  arrange(Bialko, Peptyd)

bialka_przed_ragged <- unique(trypsyna_1_przed$Bialko)
bialka_po_ragged <- unique(trypsyna_1_po$Bialko)

bialka_utracone_po_ragged <- tibble(
  Bialko = setdiff(bialka_przed_ragged, bialka_po_ragged)
) %>%
  left_join(
    fasta_df %>% select(Bialko, Nazwa_bialka),
    by = "Bialko"
  ) %>%
  arrange(Bialko)

porownanie_ragged <- tibble(
  Miara = c(
    "Peptydy unikatowe — przed filtrem ragged",
    "Peptydy unikatowe — po filtrze ragged",
    "Usunięte peptydy unikatowe",
    "Białka reprezentowane — przed filtrem ragged",
    "Białka reprezentowane — po filtrze ragged",
    "Białka utracone po filtrze ragged",
    "Peptydy panelowe — przed filtrem ragged",
    "Peptydy panelowe — po filtrze ragged",
    "Białka panelowe — przed filtrem ragged",
    "Białka panelowe — po filtrze ragged"
  ),
  Wartosc = c(
    nrow(trypsyna_1_przed),
    nrow(trypsyna_1_po),
    nrow(usuniete_ragged),
    n_distinct(trypsyna_1_przed$Bialko),
    n_distinct(trypsyna_1_po$Bialko),
    nrow(bialka_utracone_po_ragged),
    nrow(filter(trypsyna_1_przed, Bialko %in% panel_ids)),
    nrow(filter(trypsyna_1_po, Bialko %in% panel_ids)),
    n_distinct(filter(trypsyna_1_przed, Bialko %in% panel_ids)$Bialko),
    n_distinct(filter(trypsyna_1_po, Bialko %in% panel_ids)$Bialko)
  )
)

# ----------------------------------------------------------------------------
# 10. Białka panelowe bez peptydu i kontrola dłuższych kandydatów
# ----------------------------------------------------------------------------

panel_z_peptydem_po_ragged <- intersect(panel_ids, unique(trypsyna_1_po$Bialko))
panel_bez_peptydu_po_ragged <- setdiff(panel_ids, panel_z_peptydem_po_ragged)

panel_bez_peptydu_tabela <- tibble(Bialko = panel_bez_peptydu_po_ragged) %>%
  left_join(
    fasta_df %>% select(Bialko, Nazwa_bialka, Sekwencja),
    by = "Bialko"
  ) %>%
  mutate(Dlugosc_bialka = nchar(Sekwencja)) %>%
  select(Bialko, Nazwa_bialka, Dlugosc_bialka, Sekwencja)

# Globalne przypisania do 35 aa służą do oceny unikatowości dłuższych peptydów.
przypisania_trypsyna_1_do_35 <- generuj_przypisania_peptydow(
  fasta_df = fasta_df,
  enzym = "trypsyna",
  max_dpmc = 1L,
  min_len = min_dlugosc,
  max_len = max_dlugosc_kontrolna,
  usun_ragged = TRUE
) %>%
  group_by(Peptyd) %>%
  mutate(Liczba_bialek_dla_sekwencji = n_distinct(Bialko)) %>%
  ungroup()

kandydaci_dluzsi_dla_brakujacych <- przypisania_trypsyna_1_do_35 %>%
  filter(
    Bialko %in% panel_bez_peptydu_po_ragged,
    Dlugosc_peptydu > max_dlugosc,
    Dlugosc_peptydu <= max_dlugosc_kontrolna
  ) %>%
  mutate(
    Peptyd_unikatowy = Liczba_bialek_dla_sekwencji == 1L,
    Zakres_dlugosci = kategoria_dlugosci(Dlugosc_peptydu)
  ) %>%
  arrange(Bialko, Dlugosc_peptydu, Peptyd)

# Pełna diagnostyka wszystkich produktów trawienia białka bez kandydata, również tych
# odrzuconych z powodu metioniny, motywu ragged albo braku unikatowości.
wszystkie_sekwencje_do_35 <- fasta_df %>%
  mutate(
    Peptydy = map(
      Sekwencja,
      ~ trawienie_bialka(.x, enzym = "trypsyna", max_dpmc = 1L)
    )
  ) %>%
  select(Bialko, Peptydy) %>%
  unnest(Peptydy) %>%
  rename(Peptyd = Peptydy) %>%
  distinct(Bialko, Peptyd) %>%
  filter(nchar(Peptyd) <= max_dlugosc_kontrolna) %>%
  group_by(Peptyd) %>%
  mutate(Liczba_bialek_dla_sekwencji = n_distinct(Bialko)) %>%
  ungroup()

diagnostyka_brakujacych_bialek <- wszystkie_sekwencje_do_35 %>%
  filter(Bialko %in% panel_bez_peptydu_po_ragged) %>%
  left_join(
    fasta_df %>% select(Bialko, Nazwa_bialka),
    by = "Bialko"
  ) %>%
  mutate(
    Dlugosc_peptydu = nchar(Peptyd),
    Standardowe_aminokwasy = czy_standardowa_sekwencja(Peptyd),
    Zawiera_metionine = str_detect(Peptyd, "M"),
    Ragged_N = czy_ragged_N(Peptyd),
    Ragged_C = czy_ragged_C(Peptyd),
    Ragged = Ragged_N | Ragged_C,
    Peptyd_unikatowy = Liczba_bialek_dla_sekwencji == 1L,
    Powod_odrzucenia = case_when(
      Dlugosc_peptydu < min_dlugosc ~ "krótszy niż 7 aa",
      Dlugosc_peptydu > max_dlugosc_kontrolna ~ "dłuższy niż zakres kontrolny",
      !Standardowe_aminokwasy ~ "niestandardowy symbol aminokwasu",
      Zawiera_metionine ~ "zawiera metioninę",
      Ragged ~ "potencjalny ragged end",
      !Peptyd_unikatowy ~ "sekwencja współdzielona z innym białkiem",
      Dlugosc_peptydu > max_dlugosc ~ "unikatowy kandydat dłuższy niż 25 aa",
      TRUE ~ "spełnia kryteria"
    ),
    Powody_odrzucenia_wszystkie = pmap_chr(
      list(
        Dlugosc_peptydu,
        Standardowe_aminokwasy,
        Zawiera_metionine,
        Ragged,
        Peptyd_unikatowy
      ),
      function(dlugosc, standardowe, metionina, ragged, unikatowy) {
        powody <- character(0)

        if (dlugosc < min_dlugosc) {
          powody <- c(powody, "krótszy niż 7 aa")
        }
        if (dlugosc > max_dlugosc_kontrolna) {
          powody <- c(powody, "dłuższy niż zakres kontrolny")
        }
        if (!standardowe) {
          powody <- c(powody, "niestandardowy symbol aminokwasu")
        }
        if (metionina) {
          powody <- c(powody, "zawiera metioninę")
        }
        if (ragged) {
          powody <- c(powody, "potencjalny ragged end")
        }
        if (!unikatowy) {
          powody <- c(powody, "sekwencja współdzielona z innym białkiem")
        }
        if (dlugosc > max_dlugosc && dlugosc <= max_dlugosc_kontrolna) {
          powody <- c(powody, "długość przekracza podstawowy limit 25 aa")
        }

        if (length(powody) == 0) {
          "spełnia kryteria"
        } else {
          paste(unique(powody), collapse = "; ")
        }
      }
    )
  ) %>%
  arrange(Bialko, Dlugosc_peptydu, Peptyd)

# ----------------------------------------------------------------------------
# 11. Testy regresyjne i kontrole jakości
# ----------------------------------------------------------------------------

# Wartości kontrolne dla 20 663 rekordów Homo sapiens, po wyłączeniu iRT.
# Panel nie zawiera rekordu iRT, dlatego jego wartości pozostają niezmienione.
oczekiwane_przed <- tribble(
  ~Klucz,                ~Peptydy, ~Bialka, ~Panel_peptydy, ~Panel_bialka,
  "Trypsyna_0MC",       336842L,  19837L,  14610L,         739L,
  "Trypsyna_1MC_PRZED", 831849L,  20053L,  36290L,         740L,
  "LysC_0MC",           179541L,  18473L,  8439L,          724L,
  "LysC_1MC",           353613L,  18739L,  16635L,         726L
)

oczekiwane_po_ragged <- tibble(
  Peptydy = 781780L,
  Bialka = 20044L,
  Panel_peptydy = 34123L,
  Panel_bialka = 740L
)

kontrola_regresji <- map_dfr(seq_len(nrow(oczekiwane_przed)), function(i) {
  klucz <- oczekiwane_przed$Klucz[i]
  df <- unikatowe_warianty[[klucz]]
  df_panel <- df %>% filter(Bialko %in% panel_ids)

  tibble(
    Klucz = klucz,
    Miara = c("Peptydy", "Białka", "Panel_peptydy", "Panel_białka"),
    Oczekiwane = c(
      oczekiwane_przed$Peptydy[i],
      oczekiwane_przed$Bialka[i],
      oczekiwane_przed$Panel_peptydy[i],
      oczekiwane_przed$Panel_bialka[i]
    ),
    Otrzymane = c(
      nrow(df),
      n_distinct(df$Bialko),
      nrow(df_panel),
      n_distinct(df_panel$Bialko)
    )
  ) %>%
    mutate(Zgodne = Oczekiwane == Otrzymane)
})

if (!all(kontrola_regresji$Zgodne)) {
  zapisz_csv2(kontrola_regresji, "KONTROLA_REGRESJI_NIEZGODNOSC.csv")
  stop(
    "Nie odtworzono wartości kontrolnych dla zbioru Homo sapiens bez iRT. ",
    "Sprawdz plik KONTROLA_REGRESJI_NIEZGODNOSC.csv."
  )
}

# Kontrole logiczne nowego zbioru.
stopifnot(
  n_distinct(fasta_df$Bialko) == oczekiwana_liczba_rekordow_Homo_sapiens,
  n_distinct(fasta_irt$Bialko) == oczekiwana_liczba_standardow_iRT,
  !any(fasta_irt$Bialko %in% fasta_df$Bialko),
  nrow(trypsyna_1_po) == oczekiwane_po_ragged$Peptydy,
  n_distinct(trypsyna_1_po$Bialko) == oczekiwane_po_ragged$Bialka,
  nrow(filter(trypsyna_1_po, Bialko %in% panel_ids)) == oczekiwane_po_ragged$Panel_peptydy,
  n_distinct(filter(trypsyna_1_po, Bialko %in% panel_ids)$Bialko) == oczekiwane_po_ragged$Panel_bialka,
  nrow(trypsyna_1_po) <= nrow(trypsyna_1_przed),
  nrow(usuniete_ragged) == nrow(trypsyna_1_przed) - nrow(trypsyna_1_po),
  all(usuniete_ragged$Ragged),
  all(!czy_ragged(trypsyna_1_po$Peptyd)),
  all(nchar(trypsyna_1_po$Peptyd) >= min_dlugosc),
  all(nchar(trypsyna_1_po$Peptyd) <= max_dlugosc),
  all(!str_detect(trypsyna_1_po$Peptyd, "M")),
  all(czy_standardowa_sekwencja(trypsyna_1_po$Peptyd))
)

kontrola_unikatowosci_po_ragged <- przypisania_warianty[["Trypsyna_1MC"]] %>%
  filter(Peptyd %in% trypsyna_1_po$Peptyd) %>%
  count(Peptyd, name = "Liczba_przypisanych_bialek") %>%
  summarise(
    Wszystkie_maja_jedno_bialko = all(Liczba_przypisanych_bialek == 1L),
    Maksymalna_liczba_bialek = max(Liczba_przypisanych_bialek)
  )

stopifnot(kontrola_unikatowosci_po_ragged$Wszystkie_maja_jedno_bialko)

# Wpływ wyłączenia technicznego rekordu iRT na liczby FASTA.
# Wartości „przed” pochodzą z poprzedniego przebiegu, w którym iRT był liczony
# razem z białkami człowieka. Wartości panelu nie zmieniają się.
wplyw_wylaczenia_iRT <- tribble(
  ~Wariant, ~Peptydy_przed, ~Peptydy_po, ~Bialka_przed, ~Bialka_po,
  "Trypsyna, 0 DPMC", 336853L, 336842L, 19838L, 19837L,
  "Trypsyna, 1 DPMC przed filtrem ragged", 831866L, 831849L, 20054L, 20053L,
  "Trypsyna, 1 DPMC po filtrze ragged", 781797L, 781780L, 20045L, 20044L,
  "LysC, 0 DPMC", 179544L, 179541L, 18474L, 18473L,
  "LysC, 1 DPMC", 353616L, 353613L, 18740L, 18739L
) %>%
  mutate(
    Zmiana_peptydow = Peptydy_po - Peptydy_przed,
    Zmiana_bialek = Bialka_po - Bialka_przed
  )

# ----------------------------------------------------------------------------
# 12. Zapis tabel wynikowych
# ----------------------------------------------------------------------------

zapisz_csv2(kontrola_rekordow_fasta, "E1_kontrola_rekordow_FASTA.csv")
zapisz_csv2(fasta_irt, "E1_standard_iRT_WYLACZONY_Z_ANALIZY.csv")
zapisz_csv2(wplyw_wylaczenia_iRT, "E1_wplyw_wylaczenia_iRT_na_wyniki.csv")

# Końcowe podsumowanie: dokładnie cztery warianty.
zapisz_csv2(
  podsumowanie_wariantow,
  "E1_podsumowanie_4_warianty_PO_KOREKCIE.csv",
  folder_docelowy = folder_wynikowy
)
zapisz_csv2(
  podsumowanie_wariantow,
  "E1_podsumowanie_4_warianty.csv",
  folder_docelowy = folder_wynikowy
)
zapisz_csv2(
  podsumowanie_wariantow,
  "E1_podsumowanie_4_warianty_PO_KOREKCIE.csv"
)

# Osobna tabela kontrolna z Trypsyną 1 DPMC przed i po filtrze ragged.
zapisz_csv2(
  podsumowanie_wariantow_z_kontrola,
  "E1_podsumowanie_wariantow_PRZED_PO_RAGGED.csv",
  folder_docelowy = folder_wynikowy
)
zapisz_csv2(
  podsumowanie_wariantow_z_kontrola,
  "E1_podsumowanie_wariantow_PRZED_PO_RAGGED.csv"
)
zapisz_csv2(porownanie_ragged, "E1_porownanie_Trypsyna_1DPMC_przed_po_ragged.csv")
zapisz_csv2(usuniete_ragged, "E1_peptydy_usuniete_jako_ragged.csv")
zapisz_csv2(bialka_utracone_po_ragged, "E1_bialka_utracone_po_filtrze_ragged.csv")
zapisz_csv2(podsumowanie_dlugosci_fasta, "E1_rozklad_dlugosci_CALA_FASTA.csv")
zapisz_csv2(podsumowanie_dlugosci_panel, "E1_rozklad_dlugosci_PANEL_741.csv")
zapisz_csv2(dystrybucja_bialek_fasta, "E1_dystrybucja_peptydow_na_bialko_CALA_FASTA.csv")
zapisz_csv2(dystrybucja_bialek_panel, "E1_dystrybucja_peptydow_na_bialko_PANEL_741.csv")
zapisz_csv2(panel_bez_peptydu_tabela, "E1_panel_bialka_bez_peptydu_Trypsyna_1DPMC_PO_RAGGED.csv")
zapisz_csv2(kandydaci_dluzsi_dla_brakujacych, "E1_dluzsi_kandydaci_26_35aa_dla_brakujacych_bialek_panelu.csv")
zapisz_csv2(diagnostyka_brakujacych_bialek, "E1_diagnostyka_brakujacych_bialek_panelu.csv")
zapisz_csv2(kontrola_regresji, "KONTROLA_regresji_wynikow_przed_ragged.csv")
zapisz_csv2(kontrola_unikatowosci_po_ragged, "KONTROLA_unikatowosci_po_ragged.csv")

# Skoroszyt zbiorczy, latwy do przeslania i interpretacji.
write.xlsx(
  list(
    Kontrola_rekordow_FASTA = kontrola_rekordow_fasta,
    Standard_iRT_wylaczony = fasta_irt,
    Wplyw_wylaczenia_iRT = wplyw_wylaczenia_iRT,
    Podsumowanie_wariantow = podsumowanie_wariantow,
    Warianty_przed_po_ragged = podsumowanie_wariantow_z_kontrola,
    Przed_po_ragged = porownanie_ragged,
    Rozklad_dlugosci_FASTA = podsumowanie_dlugosci_fasta,
    Rozklad_dlugosci_PANEL = podsumowanie_dlugosci_panel,
    Dystrybucja_bialek_FASTA = dystrybucja_bialek_fasta,
    Dystrybucja_bialek_PANEL = dystrybucja_bialek_panel,
    Bialka_utracone_ragged = bialka_utracone_po_ragged,
    Panel_bez_peptydu = panel_bez_peptydu_tabela,
    Dluzsi_kandydaci = kandydaci_dluzsi_dla_brakujacych,
    Kontrola_regresji = kontrola_regresji
  ),
  file = file.path(folder_xlsx, "ETAP_00_01_PODSUMOWANIE_KOREKT_SKK.xlsx"),
  overwrite = TRUE
)

# Duze tabele w osobnych skoroszytach.
write.xlsx(
  usuniete_ragged,
  file = file.path(folder_xlsx, "ETAP_01_PEPTYDY_USUNIETE_RAGGED.xlsx"),
  overwrite = TRUE
)

write.xlsx(
  diagnostyka_brakujacych_bialek,
  file = file.path(folder_xlsx, "ETAP_01_DIAGNOSTYKA_BRAKUJACYCH_BIALEK_PANELU.xlsx"),
  overwrite = TRUE
)

# ----------------------------------------------------------------------------
# 13. Raport wykonania skryptu
# ----------------------------------------------------------------------------

raport <- c(
  "ETAP 00/01 — KOREKTY WG UWAG PANI PROFESOR SKK",
  paste0("Data uruchomienia: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  paste0("Plik FASTA: ", basename(plik_fasta)),
  paste0("Data FASTA wg nazwy pliku: ", data_fasta_z_nazwy_pliku),
  paste0("Źródło FASTA: ", zrodlo_fasta),
  paste0("Wszystkie rekordy w pliku FASTA: ", n_distinct(fasta_wszystkie$Bialko)),
  paste0("Rekordy Homo sapiens użyte w analizie: ", n_distinct(fasta_df$Bialko)),
  paste0("Techniczne rekordy iRT wyłączone z analizy: ", n_distinct(fasta_irt$Bialko)),
  paste0("Liczba identyfikatorow panelu: ", length(panel_ids)),
  paste0("Identyfikatory panelowe nieobecne w FASTA: ", length(panel_nieobecny_w_fasta)),
  "",
  "TRYPSYNA 1 DPMC — FILTR RAGGED",
  paste0("Peptydy unikatowe przed filtrem: ", nrow(trypsyna_1_przed)),
  paste0("Peptydy unikatowe po filtrze: ", nrow(trypsyna_1_po)),
  paste0("Usuniete peptydy unikatowe: ", nrow(usuniete_ragged)),
  paste0("Bialka reprezentowane przed filtrem: ", n_distinct(trypsyna_1_przed$Bialko)),
  paste0("Bialka reprezentowane po filtrze: ", n_distinct(trypsyna_1_po$Bialko)),
  paste0("Bialka utracone po filtrze: ", nrow(bialka_utracone_po_ragged)),
  paste0("Panel — peptydy przed filtrem: ", nrow(filter(trypsyna_1_przed, Bialko %in% panel_ids))),
  paste0("Panel — peptydy po filtrze: ", nrow(filter(trypsyna_1_po, Bialko %in% panel_ids))),
  paste0("Panel — bialka przed filtrem: ", n_distinct(filter(trypsyna_1_przed, Bialko %in% panel_ids)$Bialko)),
  paste0("Panel — bialka po filtrze: ", n_distinct(filter(trypsyna_1_po, Bialko %in% panel_ids)$Bialko)),
  "",
  paste0(
    "Bialka panelowe bez peptydu po korekcie: ",
    ifelse(
      length(panel_bez_peptydu_po_ragged) == 0,
      "brak",
      paste(panel_bez_peptydu_po_ragged, collapse = ", ")
    )
  ),
  paste0(
    "Liczba unikatowych dluzszych kandydatow 26-35 aa dla brakujacych bialek: ",
    sum(kandydaci_dluzsi_dla_brakujacych$Peptyd_unikatowy, na.rm = TRUE)
  ),
  "",
  paste0("Kontrola wartości human-only: ", ifelse(all(kontrola_regresji$Zgodne), "ZALICZONA", "NIEZALICZONA")),
  paste0("Kontrola wyłączenia iRT: ", ifelse(!any(fasta_irt$Bialko %in% fasta_df$Bialko), "ZALICZONA", "NIEZALICZONA")),
  paste0("Kontrola unikatowosci po ragged: ", ifelse(kontrola_unikatowosci_po_ragged$Wszystkie_maja_jedno_bialko, "ZALICZONA", "NIEZALICZONA")),
  "",
  paste0("Folder wynikowy: ", normalizePath(folder_wynikowy, winslash = "/", mustWork = FALSE))
)

writeLines(
  raport,
  con = file.path(folder_kontrola, "RAPORT_ETAP_00_01.txt"),
  useBytes = TRUE
)

cat("\n", paste(raport, collapse = "\n"), "\n")
cat("\nZAKONCZONO POMYSLNIE.\n")

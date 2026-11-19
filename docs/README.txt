KODY DO PRACY MAGISTERSKIEJ
Daniel Drabik
„Planowanie celowanych eksperymentów proteomicznych z wykorzystaniem środowiska R”

ZAWARTOŚĆ

01_przygotowanie_fasta_i_trawienie_in_silico.R
    Przygotowanie referencyjnego FASTA, trawienie in silico, wyznaczenie
    peptydów unikatowych oraz zastosowanie filtrów sekwencyjnych.

02_kryteria_doswiadczalne_DIA_NN.R
    Połączenie peptydów teoretycznych z raportami DIA-NN oraz zastosowanie
    kryteriów doświadczalnych. Na tym etapie powstają również listy Top 2
    i Top 3 wykorzystywane dalej.

03_listy_peptydowo_prekursorowe.R
    Uporządkowanie i eksport podstawowych list peptydowo-prekursorowych
    Top 2 i Top 3 dla całej FASTA oraz panelu 741 białek.

04_okna_retencji_i_propozycje_zamian.R
    Analiza lokalnego zatłoczenia czasu retencji oraz przygotowanie
    propozycji zamian peptydów w obrębie tego samego białka.

05_kalibracja_iRT_i_listy_rozszerzone.R
    Walidacja transferu czasu retencji przez skalę iRT oraz przygotowanie
    rozszerzonych list Top 2 i Top 3.

W pakiecie nie umieszczono osobnego skryptu służącego wyłącznie do
przygotowywania rycin do pracy.

KOLEJNOŚĆ URUCHAMIANIA

Skrypty należy uruchamiać kolejno od 01 do 05. Każdy etap tworzy własny
katalog wynikowy oznaczony datą i czasem uruchomienia. Kolejne skrypty
odnajdują wyniki wcześniejszych etapów w bieżącym katalogu roboczym
i jego podkatalogach.

Najwygodniej utworzyć osobny katalog projektu, umieścić w nim pięć skryptów
oraz wszystkie zewnętrzne pliki wejściowe, ustawić ten katalog jako katalog
roboczy w RStudio i uruchamiać skrypty po kolei.

PAKIETY R

Do wykonania całego pipeline'u wykorzystywane są:
tidyverse, readxl, openxlsx, ggplot2, arrow, scales.

ZEWNĘTRZNE PLIKI WEJŚCIOWE

1. Referencyjny FASTA Homo sapiens z dołączonym rekordem standardów iRT:
   2025_11_19_homo_sapiens_20663_UP000005640_9606_iRT.fasta

   Skrypt 01 akceptuje również nazwę:
   2025_11_19_homo_sapiens_20663_UP000005640_9606_iRT(1).fasta

2. Panel 741 białek:
   ALL_candidates_741.xlsx

   Skrypt 01 akceptuje również:
   ALL_candidates_741_REKONSTRUKCJA_ID.xlsx

3. Cztery główne raporty DIA-NN wygenerowane z opcją --export-quant:
   report_2025_30_min_export_quant_added.parquet
   report_2026_30_min_export_quant_added.parquet
   report_2025_60_min_r1_export_quant_added.parquet
   report_2025_60_min_r2_export_quant_added.parquet

Pliki RAW i mzML nie są bezpośrednim wejściem do tych skryptów R. Zostały
wykorzystane na wcześniejszym etapie przetwarzania danych w DIA-NN. Do
ponownego uruchomienia przedstawionego pipeline'u R potrzebne są wyłącznie
FASTA, panel 741 białek i cztery wymienione raporty DIA-NN.

ZALEŻNOŚCI POMIĘDZY ETAPAMI

Etap 01:
    wejście zewnętrzne -> FASTA + panel 741 białek

Etap 02:
    wejście -> wyniki etapu 01 + panel 741 białek + 4 raporty DIA-NN

Etap 03:
    wejście -> wyniki etapu 02

Etap 04:
    wejście -> wyniki etapu 03 + wybrane tabele z etapu 02

Etap 05:
    wejście -> wyniki etapu 03 + wybrane tabele i dane pośrednie z etapu 02

UWAGA

W kodach pozostawiono automatyczne wyszukiwanie wyników wcześniejszych etapów,
ale ograniczono je do katalogu roboczego i katalogu zawierającego skrypty.
Dzięki temu skrypty nie wybierają przypadkowych kopii plików z katalogów
Downloads lub Desktop.

Wersja załącznikowa zachowuje wykonywalne instrukcje R z wersji użytych do
końcowych obliczeń. Zmiany dotyczą komentarzy oraz nazw samych plików .R;
logika analityczna i działające mechanizmy wyszukiwania danych nie zostały
zmienione.

# I. ЧАСТЬ. Получение идентификаторов интересующих генов, их транскриптов и геномных диапазонов GRanges

## 1.1. Загрузка необходимых библиотек
## 1.1.1. Установка библиотек, если ранее не были установлены
packages_bioconductor <- c("org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene", "BSgenome.Hsapiens.UCSC.hg38", "GenomicRanges", "VariantAnnotation")
for(package in packages_bioconductor){
  if (!requireNamespace(package, quietly = TRUE)) {
    BiocManager::install(package)
  }
}

packages_cran <- c("dplyr", "httr", "jsonlite", "tidyr")
for(package in packages_cran){
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
}

## 1.1.2. Подключение библиотек, используемые в ходе работы:
library(dplyr)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(VariantAnnotation)

## 1.2. Получение идентификаторов генов по их символам и создание объекта GRanges для последующей фильтрации данных:
### 1.2.1. Получение идентификаторов генов Entrez ID по их символу из базы "org.Hs.eg.db"
genes_symbol <- c("SLC3A1", "TTN", "GUF1", "TMEM70", "CTU2", "MKKS")
genes_id <- AnnotationDbi::select(org.Hs.eg.db, keys=genes_symbol, keytype="SYMBOL", columns=c("ENTREZID", "ENSEMBL", "SYMBOL"))
genes_ids <- genes_id$ENTREZID

### 1.2.2. Получение списка транскриптов генов по идентификаторам генов Entrez ID из базы "TxDb.Hsapiens.UCSC.hg38.knownGene"
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_txGRL <- transcriptsBy(txdb, "gene")[genes_ids]

### 1.2.3. Формирование таблицы (DataFrame) из списка генов и транскриптов
genes_tx_list <- lapply(names(genes_txGRL), function(id) {
  gr <- genes_txGRL[[id]]
  df <- as.data.frame(gr)
  df$GENE_SYMBOL <- genes_id$SYMBOL[match(id, genes_id$ENTREZID)]
  return(df)
})
genes_txDF <- do.call(rbind, genes_tx_list)

### 1.2.4. Создание объекта GRanges по списку транскриптов
filter_by_tx <- list(tx_name=genes_txDF$tx_name)
genes_txGR <- transcripts(txdb, columns=c("tx_id", "tx_name"), filter=filter_by_tx, use.names=TRUE)

## 1.3. Импорт данных файла .vcf и фильтрация вариантов по найденным транскриптам(геномным координатам)
### 1.3.1. Импорт файла .vcf в рабочую область R
vcf_data <- readVcf("~/R/bind/data/fin_filt_FDP_PID25001.vcf", genome="hg38")

### 1.3.2. Фильтрация vcf-данных по геномным диапазонам объектов GRanges
vcf_by_tx <- subsetByOverlaps(vcf_data, genes_txGR)

### 1.3.3. Преобразование результата в DataFrame
vcf_by_txDF <- as.data.frame(rowRanges(vcf_by_tx))

## 1.4. Анализ генотипов и глубин прочтений. Извлечение и добавление генотипа и глубину прочтений в итоговую таблицу.
### 1.4.1. Получение данных о генотипах (GT) и глубине прочтений (DP)
matGT <- geno(vcf_by_tx)$GT
matDP <- geno(vcf_by_tx)$DP
num_samples <- ncol(matGT)

### 1.4.2. Добавление колонки с генотипами и описанием генотипов для каждого образца
for (i in 1:num_samples) {
  sample_name <- colnames(matGT)[i]  # Имя образца
  # Добавление числового обозначения генотипа
  vcf_by_txDF[[paste(sample_name, "_GT", sep = "")]] <- matGT[, i]
  # Добавление словесного описания генотипа на русском языке
  description <- case_when(
    matGT[, i] == "0/0" ~ "Гомозиготный дикий тип",
    matGT[, i] == "0/1" | matGT[, i] == "1/0" ~ "Гетерозиготный",
    matGT[, i] == "1/1" ~ "Гомозиготный мутантный тип",
    matGT[, i] == "1/2" ~ "Полимерический генотип",
    TRUE ~ "Неопределённое состояние"
  )
  vcf_by_txDF[[paste(sample_name, "_GT_Description", sep = "")]] <- description
  # Добавление числового обозначения глубины прочтения
  vcf_by_txDF[[paste(sample_name, "_DP", sep = "")]] <- matDP[, i]
}

## 1.5. Связывание мутаций с соответствующими ей транскриптами и генами
### 1.5.1. Добавление колонки с именем транскрипта и символом гена
vcf_by_txDF$tx_name <- NA
vcf_by_txDF$GENE_SYMBOL <- NA

### 1.5.2. Прохождение по каждой мутации и назначение транскрипта и гена
for (i in 1:nrow(vcf_by_txDF)) {
  # Координата мутации
  mut_pos <- vcf_by_txDF[i, "start"]
  
  # Проверка, какая мутация попадает в диапазон транскрипта
  matching_tx <- genes_txDF[(genes_txDF$start <= mut_pos & genes_txDF$end >= mut_pos), ]
  
  # Если есть нахождение подходящего транскрипта, фиксируется его имя и символ гена
  if (nrow(matching_tx) > 0) {
    vcf_by_txDF[i, "tx_id"] <- matching_tx$tx_id[1]
    vcf_by_txDF[i, "tx_name"] <- matching_tx$tx_name[1]
    vcf_by_txDF[i, "GENE_SYMBOL"] <- matching_tx$GENE_SYMBOL[1]
  }
}

### 1.5.3. Вывод результата
print(vcf_by_txDF)

## Следуя данному руководству, вы получите полную картину ваших генетических данных, сформировав удобную таблицу с необходимой информацией о мутациях и соответствующих генах.

# II. ЧАСТЬ. Использование базы данных gnomAD
# gnomAD (Genome Aggregation Database) — крупная база данных, содержащая миллионы генетических вариантов среди населения Европы и других регионов мира. Этот раздел посвящен извлечению данных по целевым генам и вычислению частот встречаемости аллелей в европейской популяции (без учета финнов, сокращенно NFE).

## 2.1. Подключение библиотек для запросов и обработки данных и подготовка окружения
library(httr)        # HTTP запросы
library(jsonlite)    # Парсер JSON
library(dplyr)       # Манипуляции с данными
library(tidyr)       # Обработка сложных структур данных


## 2.2. Формирование запроса к API gnomAD.  Определение функций для получения данных по генам

### 2.2.1. Создание функции для получения списка вариантов для указанного гена
fetch_gene_variants <- function(gene_symbol) {
  # Графический запрос (GraphQL) для извлечения информации о вариантах
  query <- paste0('query VariantsInGene {',
                  '  gene(gene_symbol: "', gene_symbol, '", reference_genome: GRCh38) {',
                  '    variants(dataset: gnomad_r4) {',
                  '      variant_id, chrom, pos, ref, alt, rsids, transcript_id, transcript_version, hgvs, hgvsc, hgvsp, consequence, flags, ',
                  '      joint { ac, an, populations { id, ac, an, homozygote_count, hemizygote_count } }',
                  '    }',
                  '  }',
                  '}')
  
  # Параметры для запроса
  url <- "https://gnomad.broadinstitute.org/api/"
  headers <- add_headers(Content_Type = "application/json")
  body <- list(query = query)
  
  # Отправка запроса
  response <- httr::POST(url, body = body, encode = "json", .headers = headers)
  
  # Проверка статуса ответа
  if (response$status_code != 200) {
    stop("Ошибка при отправке запроса:", content(response, "text"))
  }
  
  # Чтение и обработка данных
  raw_data <- content(response, "text")
  result <- fromJSON(raw_data)
  df <- as.data.frame(result$data$gene$variants, stringsAsFactors = FALSE)
  df$gene_symbol <- gene_symbol
  
  return(df)
}

## 2.3. Сбор данных по выбранным генам
### 2.3.1. Запрос данных для генов из списка
gene_symbols <- c("SLC3A1", "TTN", "GUF1", "TMEM70", "CTU2", "MKKS")
gene_var_gad <- lapply(gene_symbols, fetch_gene_variants)

## 2.4. Получение единой таблицы по интересующей информации
### 2.4.1. Объединение результатов в единую таблицу
gene_var_gadDF <- dplyr::bind_rows(gene_var_gad)

### 2.4.2. Расчёт частоты аллеля (AF) на основе общего числа аллелей
gene_var_gadDF$joint_af <- with(gene_var_gadDF, ifelse(!is.na(joint$ac) & !is.na(joint$an) & joint$an > 0, joint$ac / joint$an, NA))

### 2.4.3. Получение таблицы по популяции "nfe" (европейские нефинны)
nfe_pop_data <- lapply(seq_along(gene_var_gadDF$joint$populations), function(i) {
  pops <- gene_var_gadDF$joint$populations[[i]]
  filtered_pops <- pops %>% filter(id == "nfe") %>% dplyr::select(ac_nfe = ac, an_nfe = an, homozygote_count_nfe = homozygote_count, hemizygote_count_nfe = hemizygote_count)
  filtered_pops$variant_id <- gene_var_gadDF$variant_id[i]
  filtered_pops
}) %>% bind_rows()

### 2.4.4. Объединение данных по популяции с общей таблицей
gene_var_nfe_gadDF <- left_join(gene_var_gadDF, nfe_pop_data, by = "variant_id")

### 2.4.5. Выделение интересующих колонок из итоговой таблицы
gene_var_nfe_gadDF_sub <- subset(gene_var_nfe_gadDF, select = c(variant_id, chrom, pos, ref, alt, rsids, transcript_id, transcript_version, hgvs, consequence, gene_symbol, joint_af, ac_nfe, an_nfe, homozygote_count_nfe, hemizygote_count_nfe))

### 2.4.6. Присвоение имен строкам в таблице
vcf_by_txDF$variant <- rownames(vcf_by_txDF)

## В итоге получаем подробную таблицу по каждому варианту, включающую геномные координаты, описание последствий мутации, частоту аллелей и популяционную статистику по выбранной группе (европейской популяции без финнов).
## Теперь у вас есть готовый набор данных для дальнейшего анализа, включая частоты аллелей и детальную информацию о влиянии выявленных генетических вариантов.

# III. ЧАСТЬ. Стандартизация данных и удаление дублей
# Цель данного этапа — привести данные к единому виду и устранить повторяющиеся записи, возникающие при слиянии двух источников.

## 3.1. Приведение идентификаторов вариантов к формату gnomAD
### 3.1.1. Создание функции для преобразования идентификаторов генетических вариантов в формат gnomAD
format_variant <- function(variant_string) {
  # Парсим входящую строку
  parts <- strsplit(variant_string, ":|_", fixed = FALSE)[[1]]
  chrom <- gsub("chr", "", parts[1])   # Убираем префикс "chr"
  pos <- parts[2]                      # Геномная позиция
  alleles <- strsplit(parts[3], "/")[[1]]  # Альтернативные нуклеотиды
  ref <- alleles[1]                   # Реверенс-нуклеотид
  alt <- alleles[2]                   # Альтернативный нуклеотид
  # Составляем итоговую строку в формате gnomAD
  paste(chrom, pos, ref, alt, sep = "-")
}

### 3.1.2. Применение функции к таблице для формирования колонки "variant_id"
vcf_by_txDF$variant_id <- sapply(vcf_by_txDF$variant, format_variant)

## 3.2. Соединение таблиц по идентификаторам вариантов. Объединение основной таблицы с результатами анализа по соответствию идентификаторов вариантов.
### 3.2.1. Левостороннее соединение таблиц по ключу "variant_id"
vcf_by_tx_gadDF <- vcf_by_txDF %>% left_join(gene_var_nfe_gadDF_sub, by = "variant_id")

## 3.3. Оптимизация структуры данных и устранение избыточности.
### 3.3.1. Преобразование столбца "transcript_id" путём объединения с номером версии транскрипта
vcf_by_tx_gadDF$transcript_id <- paste(vcf_by_tx_gadDF$transcript_id, vcf_by_tx_gadDF$transcript_version, sep = ".")

### 3.3.2. Удаление ненужных столбцов и повторяющихся столбцов
vcf_by_tx_gadDF <- vcf_by_tx_gadDF[, -which(names(vcf_by_tx_gadDF) %in% c("ref", "alt", "chrom", "pos", "transcript_version", "variant_id", "gene_symbol"))]

## 3.4. Вычисление частоты аллелей в целевой популяции. Расчет относительной частоты аллелей (allele frequency) в евразийской популяции (кроме финнов)
### 3.4.1. Добавление новой колонки с частотой аллелей (AF) для популяции "nfe"
vcf_by_tx_gadDF$af_nfe <- vcf_by_tx_gadDF$ac_nfe / vcf_by_tx_gadDF$an_nfe

## Итоговая таблица включает исчерпывающие данные по каждому варианту с указанием координат, типов мутаций, генотипов и статистики по исследованной популяции.

# IV. ЧАСТЬ. Фильтрация вариантов по характеристикам
# На данном этапе проводится детализированная фильтрация вариантов на основании частоты аллелей, типа мутации и гетерозиготности. Цель — оставить наиболее важные и потенциально значимые генетические события.

## 4.1. Фильтрация данных по специфичным условиям,учитывая наличие гетерозиготных состояний, низкую частоту аллелей и типы значимых эффектов мутаций
### 4.1.1. Фильтрация данных по условиям и выделение важных полей
vcf_by_tx_gad_filtDF <- vcf_by_tx_gadDF %>%
  # Отбор только гетерозиготных вариантов с низкой частотой (< 1%)
  filter(SampleName_GT_Description == "Гетерозиготный", af_nfe < 0.01) %>%
  # Выбор значимых типов мутаций
  filter(consequence %in% c("missense_variant", "stop_gained", "frameshift_variant", "start_lost")) %>%
  # Оставляем лишь критически важные столбцы
  dplyr::select(
    GENE_SYMBOL, tx_name, variant, QUAL, GT=SampleName_GT, DP=SampleName_DP, hgvs, rsids,
    GT_Description=SampleName_GT_Description, af_nfe, consequence
  )

## 4.2. Предварительный просмотр итоговых данных
### 4.2.1. Просмотр первую часть полученной фильтрованной таблицы
head(vcf_by_tx_gad_filtDF)

### 4.2.2. Просмотр полностью всей таблицы
print(vcf_by_tx_gad_filtDF)

## Итоговая таблица содержит отобранные варианты, удовлетворяющие заданным критериям отбора, что позволяет сосредоточиться на ключевых генетических событиях, подлежащих дальнейшему изучению.

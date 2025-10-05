# bio_utils/example_data.py
# Структура словаря: {имя_последовательности: (нуклеотиды, строка_качества)}
FASTQ_SEQS = {
    '@SEQ_ID_1': (
                  'GGTCGAAAACTGATCGGCTGT',
                  '!"#$%&\'()*+,-./0123456'
                 ),
    # Среднее качество, средний GC
    '@SEQ_ID_2': (
                  'TCTCCATGGCTTCGCTCTGCA',
                  'ABCDEFGHIJJKLMNOPQRST'
                 ),
    # Высокое качество, высокий GC
    '@SEQ_ID_3': (
                  'AGCTTGAAACGT',
                  '@@@@@@@@@@@@'
                 ),
    # Идеальное качество, средний GC
    '@SEQ_ID_4': (
                  'GATTACA',
                  '1234567'
                 ),
    # Короткая, среднее качество
    '@SEQ_ID_5': (
                  'TTTTTTTTTTTTTTTTTTTTT',
                  'FFFFFFFFFFFFFFFFFFFFF'
                 ),
    # Низкий GC, высокое качество
}

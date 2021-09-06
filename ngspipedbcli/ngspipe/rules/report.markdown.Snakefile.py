rule report_html_by_snakemake:
    input:
    output: join("../report", "report.html")
    params:
        rmd = "index.Rmd"
    script:
        "/Users/zhangxuan/Work/Projects/bioinformatics/RNAPhoenix11.14/report/index.Rmd"
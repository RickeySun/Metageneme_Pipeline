# -*- encoding: utf-8 -*-
"""
宏基因组测序fastq质检
使用前提：所有的原始fastq.gz文件需要在当前工作路径下(./)
使用方法：直接运行python fastqcResult.py
结果：包含每个样本的独立质检结果（例如：03A0136_3.raw_2_fastqc.html），以及汇总结果（multiqc_report.html）。/
    同时生工测序样本结果需满足 1、数据量2G。2、GC含量小于65%。汇总结果在fastqc_summary.csv中。
    如质量不满足需反馈生工，重新测序。
"""

import os
import subprocess
import pandas as pd
from bs4 import BeautifulSoup
from multiprocessing import Pool


def run_command(command):
    print(f"Running command: {command}")
    process = subprocess.Popen(
        command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    # 实时打印输出内容
    while True:
        output = process.stdout.readline()
        error = process.stderr.readline()

        if output:
            print(output.strip())
        if error:
            print(error.strip())

        if output == "" and error == "" and process.poll() is not None:
            break

    process.wait()  # 等待命令完成
    if process.returncode != 0:
        print(f"Error: Command failed with return code {process.returncode}")
    else:
        print("Command completed successfully!")

# 单个样本运行 FastQC
def run_fastqc(filename, input_dir, threads):
    html_filename = filename.replace('.fastq.gz', '_fastqc.html')
    if html_filename not in os.listdir(input_dir):
        run_command(f"fastqc {os.path.join(input_dir, filename)} -t {threads}")
    else:
        print(f"FastQC HTML file for {filename} already exists. Skipping FastQC step.")

# 运行 FastQC 和 MultiQC
def run_fastqc_and_multiqc(input_dir, threads=10):
    fastq_files = [f for f in os.listdir(input_dir) if f.endswith(".fastq.gz")]

    # 使用多进程并行运行 FastQC
    with Pool(processes=threads) as pool:
        pool.starmap(run_fastqc, [(f, input_dir, 1) for f in fastq_files])

    # 运行 MultiQC
    run_command(f"multiqc {input_dir}")

# 解析 FastQC HTML 文件提取信息
def parse_fastqc_html(input_dir):
    results = []

    for filename in os.listdir(input_dir):
        if filename.endswith("_fastqc.html"):
            filepath = os.path.join(input_dir, filename)

            # 打开并解析HTML文件
            with open(filepath, "r", encoding="utf-8") as file:
                soup = BeautifulSoup(file, "html.parser")

                # 提取Filename
                filename_tag = soup.find("td", string="Filename")
                filename_value = filename_tag.find_next("td").text if filename_tag else "N/A"

                # 提取Total Bases
                bases_tag = soup.find("td", string="Total Bases")
                bases_value = bases_tag.find_next("td").text if bases_tag else "N/A"

                # 提取%GC
                gc_tag = soup.find("td", string="%GC")
                gc_value = gc_tag.find_next("td").text if gc_tag else "N/A"

                # 将结果存入列表
                results.append({
                    "Filename": filename_value,
                    "Total Bases": bases_value,
                    "%GC": gc_value
                })

    # 返回提取的结果
    return pd.DataFrame(results)


def main():
    input_dir = "./"  # 当前目录
    output_file = "fastqc_summary.csv"

    # 运行 FastQC 和 MultiQC
    print("Running FastQC and MultiQC...")
    run_fastqc_and_multiqc(input_dir)

    # 解析 HTML 文件
    print("Parsing FastQC HTML files...")
    summary_df = parse_fastqc_html(input_dir)

    # 保存结果到 CSV 文件
    summary_df.to_csv(output_file, index=False)
    print(f"Summary saved to {output_file}")


if __name__ == "__main__":
    main()

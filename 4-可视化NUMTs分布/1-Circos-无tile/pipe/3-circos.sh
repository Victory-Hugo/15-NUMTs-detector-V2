#!/bin/bash
set -euo pipefail
#* 需要先安装Circos工具,可以使用conda或者下载安装
#* conda install -c conda-forge libjpeg-turbo perl-gd
#* conda install -c bioconda circos

# 按脚本自身位置定位项目目录，避免硬编码路径失效。
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
CONF_DIR="${PROJECT_DIR}/conf"
CONF_FILE="circos_allNUMTs.conf"
CIRCOS_BIN="${CIRCOS_BIN:-}"

if [[ ! -d "${CONF_DIR}" ]]; then
    echo "错误: 配置目录不存在: ${CONF_DIR}" >&2
    exit 1
fi

if [[ ! -f "${CONF_DIR}/${CONF_FILE}" ]]; then
    echo "错误: 配置文件不存在: ${CONF_DIR}/${CONF_FILE}" >&2
    exit 1
fi

if [[ -z "${CIRCOS_BIN}" ]]; then
    if command -v circos >/dev/null 2>&1; then
        CIRCOS_BIN="$(command -v circos)"
    elif [[ -n "${CONDA_PREFIX:-}" && -x "${CONDA_PREFIX}/bin/circos" ]]; then
        CIRCOS_BIN="${CONDA_PREFIX}/bin/circos"
    else
        echo "错误: 未找到 circos 可执行文件。请先激活包含 circos 的 Conda 环境，或通过 CIRCOS_BIN 指定路径。" >&2
        exit 1
    fi
fi

MODULE_CHECK_LOG="$(mktemp)"
trap 'rm -f "${MODULE_CHECK_LOG}"' EXIT

"${CIRCOS_BIN}" -modules >"${MODULE_CHECK_LOG}" 2>&1 || true
if grep -q '^missing' "${MODULE_CHECK_LOG}"; then
    echo "错误: 当前 Circos 运行环境缺少 Perl 依赖模块，无法正常绘图。" >&2
    echo "缺失模块如下:" >&2
    grep '^missing' "${MODULE_CHECK_LOG}" >&2
    echo "请先在当前环境中补齐依赖后再运行。" >&2
    exit 1
fi

#! 使用相对路径，因为circos配置太难了
cd "${CONF_DIR}"

# 执行前配置好conf文件
"${CIRCOS_BIN}" -conf "${CONF_FILE}"

# 绿色输出：成功提示
echo -e "\e[32m===============================================\e[0m"
echo -e "\e[32mCircos可视化NUMTs分布完成\e[0m"
echo -e "\e[32m输出文件在output目录下\e[0m"
echo -e "\e[32m===============================================\e[0m"

# 青色输出：配置修改提示
echo -e "\e[36m如果需要修改图片格式,请修改conf文件中的image部分\e[0m"
echo -e "\e[36m如果需要修改链接样式,请修改conf文件中的links部分\e[0m"
echo -e "\e[36m如果需要修改染色体样式,请修改conf文件中的ideogram部分\e[0m"
echo -e "\e[36m如果需要修改颜色样式,请修改conf文件中的custom_colors.conf部分的llt_color部分，使用RGB编码\e[0m"
echo -e "\e[36m如果需要修改链接规则,请修改conf文件中的rules部分\e[0m"
echo -e "\e[36m如果需要修改其他参数,请参考Circos官方文档\e[0m"
echo -e "\e[36mCircos官方文档: http://circos.ca/documentation/\e[0m"

# 黄色输出：输入输出文件修改提示
echo -e "\e[33m===============================================\e[0m"
echo -e "\e[33m如果需要修改输入文件,请修改 conf 文件中的 file = ../input/circos_SLE.txt\e[0m"
echo -e "\e[33m如果需要修改输出文件,请修改 file* = ../output/CircosPlot_allNUMT\e[0m"
echo -e "\e[33m===============================================\e[0m"

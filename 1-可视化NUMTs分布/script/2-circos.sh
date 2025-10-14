#!/bin/bash
#* 需要先安装Circos工具,可以使用conda或者下载安装
#* conda install -c bioconda circos

# 切换到配置文件目录（因为配置文件中使用相对路径）
SCRIPT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/1-可视化NUMTs分布"

#! 使用相对路径，因为circos配置太难了
cd ${SCRIPT_DIR}/conf/

# 执行前配置好conf文件
circos -conf circos_allNUMTs.conf

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
echo -e "\e[33m如果需要修改输入文件,请修改 file = ../output/circos.txt\e[0m"
echo -e "\e[33m如果需要修改输出文件,请修改 file* = ../output/CircosPlot_allNUMT\e[0m"
echo -e "\e[33m===============================================\e[0m"

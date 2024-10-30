# write awk script to extract Na begin row

# 2024-10-21 zss: 看上去只能提取Na开头的行，但是这样提取出来并不是dump文件！
cat > extract_dump.awk << EOF
BEGIN {
    FS = " "
    OFS = " "
    print "Na"
}
/^Na/ {
    print \$0
}
EOF

awk -e extract_dump.awk $1 > $2
# write awk script to extract Na begin row
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
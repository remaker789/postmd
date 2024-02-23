# after build the docs, run this powershell script.
Remove-Item -Path "html/*" -Recurse
Remove-Item -Path "../postmd-manual/*" -Recurse

sphinx-multiversion.exe docs html

# 确保目标文件夹存在，如果不存在则创建它
New-Item -ItemType Directory -Path "../postmd-manual/" -Force -ErrorAction SilentlyContinue

# 移动源目录及其所有子目录和文件到目标目录
Move-Item -Path "html/*" -Destination "../postmd-manual/"


cd ../postmd-manual/
New-Item -ItemType File -Path ".nojekyll"

# 获取当前时间并格式化
$timestamp = Get-Date -Format "yyyy-MM-dd HH:mm"

git add .
git commit -m "update at $timestamp"
git push -u origin main


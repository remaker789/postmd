from PostMD.FileOperation import *
test = Dir(r"D:\OneDrive\Github\MyRepos\PostMD")
print(test)
out = test.get_path(pattern="py",mode="suffix",recursive_hier='1')
print(out)
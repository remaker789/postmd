## 2022-09-21

1. mpl的latex应该是只能用Times。如果以后需要改成非衬线字体，再一个个改吧。。。
2. 对于不适用latex进行渲染，而是使用mathtext进行渲染，如果需要自定义字体，一方面需要使用 `mathtext.fontset: custom`, 然后由于自定义字体的数学格式（如上下标间距）没有很好地定义，因此需要进行修改。对于新版本的mpl，需要进行以下修改：

   ```python
   from matplotlib.mathtext import _mathtext as mathtext
   mathtext.FontConstantsBase = mathtext.ComputerModernFontConstants
   ```

   修改之后上下标间距就正常多了。

   对于Type42，输出pdf图片时，会出现**warning："meta NOT subset; don't know how to subset; dropped"**， 目前还不知道什么问题。

***建议无衬线字体使用Type42-Arial-totally.mplstyle；衬线字体使用Type42-Times-totally.mplstyle！***

**目前mplstyle已经在matplotlib-3.6.0中测试成功！**

## 2023-08-15

适配matplotlib 3.7.2





## 2024-12-09

从scienceplot复制了styles

# class Dir:
#     def __init__(self, dir):
#         self.dir = dir
#         self.Path = Path(self.dir)
#         if not self.Path.is_dir():
#             raise ValueError("The keyword dir must be a directory")
    
#     # @staticmethod             
#     # def dig_next_depth(Path, pattern):
#     #     pattern = '*/'+pattern # bug: 不知道*/对应linux系统是否适用。
#     #     return Path.glob(pattern) # 返回dir下面所有files
                
        
#     def get_path(self, 
#                  pattern:str = None, 
#                  mode = 'full', 
#                  rex=False, 
#                  recursive_hier='all',
#                  out = "str") -> list :
#         """get filepaths according to pattern

#         Args:
#             pattern (str, optional): pattern used to math the filename in Dir. Defaults to None.
#             mode (str, optional): the match criterion of pattern. Defaults to 'full'. 
#                                   full match if mode=='all';
#                                   match prefix if mode=='prefix';
#                                   match suffix if mode=="suffix';
#                                   fnmatch mode if mode=="custom". see https://docs.python.org/zh-cn/3/library/fnmatch.html#module-fnmatch
#             re (bool, optional): use RegularExpression if re is True, else fnmatch mode. Defaults to False.
#             recursive_hier (str, optional): The recursive hierarchy to dig. Defaults to 'all'. 
#                                             all depth will be searched if recursive_hier='all';
#                                             the specified depth will be searched if recusive_hier is int (the base depth is 0);
                                            
                                            
#             out (str, optional): _description_. Defaults to "str".

#         Raises:
#             ValueError: _description_
#             ValueError: _description_
#             ValueError: _description_

#         Returns:
#             list: _description_
#         """        
#         if pattern is None:
#             raise ValueError("A pattern is desired to get the filepaths")
        
#         # fnmatch if re is False; 
#         if not rex:
#             if mode == "full":       # full match
#                 pattern = pattern
#             elif mode == "prefix":   # match prefix
#                 pattern = pattern+"*"
#             elif mode == "suffix":   # match suffix
#                 pattern = '*'+pattern    
#             elif mode == "custom":   # custom pattern
#                 pattern = pattern
#             else:
#                 raise ValueError("input mode is not supported!")
#         else: # re的还没有做。。。
#             mode == None
#             print("RegularExpression is used!")
            
#         self.short_pattern = {"pattern": pattern, "re": rex}
#         self.re = rex
        
#         if recursive_hier == 'all':
#             pattern = '**/' + pattern # 递归所有目录
#             _paths = self.Path.glob(pattern) # 返回dir下面所有files

#         elif isinstance(recursive_hier, int):
#             pattern = recursive_hier*'*/'+pattern # bug: 不知道*/对应linux系统是否适用。
#             _paths = self.Path.glob(pattern) # 返回dir下面所有files

#         elif isinstance(recursive_hier, str): # 目前最大允许递归到第9层目录
#             _paths = []
#             if re.match('\d\*\d',recursive_hier): # 从first_depth检索到end_depth
#                 first_depth = int(recursive_hier[0])
#                 end_depth = int(recursive_hier[-1])
#                 for depth in range(first_depth, end_depth+1):
#                     depth_pattern = depth*'*/'+pattern
#                     _paths += self.Path.glob(depth_pattern)
                    
#             elif re.match('\*\d', recursive_hier): # 从第0层检索到end_depth
#                 first_depth = 0
#                 end_depth = int(recursive_hier[-1])
#                 for depth in range(first_depth, end_depth+1):
#                     depth_pattern = depth*'*/'+pattern
#                     _paths += self.Path.glob(depth_pattern)
#             elif re.match('\d\*',recursive_hier): # 从first_depth检索到最后
#                 first_depth = int(recursive_hier[0])
#                 end_depth = 100 # 将最大层设置成100，只要没超过100层，应该就可以吧。。。
#                 for depth in range(first_depth, end_depth):
#                     depth_pattern = depth*'*/'+pattern
#                     _paths += self.Path.glob(depth_pattern)
#             elif re.match('\d+', recursive_hier): # 只检索第n层
#                 depth_pattern = int(recursive_hier)*'*/'+pattern
#                 _paths += self.Path.glob(depth_pattern) # 返回dir下面所有files
#             else:
#                 raise ValueError("the fommat of input string is not supported")
            
#         # elif 如果recursive_hier是可迭代对象，对里面每一个数值进行迭代
        
#         else:
#             raise ValueError("Input recursive_hier is not supported.")
        
        
        
        
#         # 选出files
#         if out == 'str':
#             files = [str(_path) for _path in _paths if _path.is_file() ]
#             # dirs = [str(_path) for _path in _paths if _path.is_dir() ]
#         else: 
#             files = [_path for _path in _paths if _path.is_file() ]
#             # dirs = [_path for _path in _paths if _path.is_dir() ]
#         return files
    
# coding: utf-8

# We will use semantic version (major, minor, patch)
__main_version_tuple__ = main_version_tuple = (0, 0, 0)
__pre_version_tuple__ = pre_version_tuple = ('alpha',1) # or None or an empty tuple, the first char is alpha, beta, rc, etc...
if __pre_version_tuple__:
    __version_tuple__ = version_tuple = tuple(list(__main_version_tuple__) + list(__pre_version_tuple__))
    __version__ = version = '.'.join([str(i) for i in __main_version_tuple__]) + f'-{__pre_version_tuple__[0]}' + '.'.join([str(i) for i in __pre_version_tuple__[1:]])
else:
     __version__ = version = '.'.join([str(i) for i in __main_version_tuple__])
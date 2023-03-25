import os
import inspect

path = os.path.realpath(__file__)
src_path = path.rsplit('/', 1)[0]
print(src_path)


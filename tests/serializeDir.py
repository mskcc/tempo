#!/usr/bin/env python3
"""
Module for serializing directory file lists for use with testing
"""
import os
import sys
import hashlib
import json
from pathlib import Path

class DirSerializer(object):
    """
    Serializes the listing of all files in a directory tree
    """
    def __init__(self, root, exclude_dirs = ['work'], exclude_exts = ['.log']):
        self.root = root # path to start from
        self.exclude_dirs = exclude_dirs

        self.data = {}

        for dirname, dirs, files in os.walk(root):
            # remove excluded dirs and files from the list
            dirs[:] = [ d for d in dirs if d not in self.exclude_dirs ]
            files[:] = [ fn for fn in files if any(not fn.endswith(ext) for ext in exclude_exts) ]

            for file in files:
                path = os.path.join(dirname, file)
                key = str(Path(path).relative_to(root)) # dont include the input dir name in the key
                self.data[key] = {}
                self.data[key]['md5'] = self.md5sum(path)
                self.data[key]['size'] = os.path.getsize(path)
                if self.isText(path):
                    self.data[key]['lines'] = self.numLines(path)

    def json(self, **kwargs):
        """
        print(DirSerializer('/paths').json(indent = 4))
        """
        return(json.dumps(self.data, **kwargs))

    def md5sum(self, filepath, chunk_num_blocks=128):
        file_hash=hashlib.md5()
        with open(filepath,'rb') as f:
            for chunk in iter(lambda: f.read(chunk_num_blocks * file_hash.block_size), b''):
                file_hash.update(chunk)
        return(file_hash.hexdigest())

    def numLines(self, filepath):
        count = 0
        with open(filepath) as f:
            for i, _ in enumerate(f):
                count += 1
        return(count)

    def isText(self, filepath):
        try:
            with open(filepath, "r") as f:
                for line in f:
                    break
            return(True)
        except UnicodeDecodeError:
            return(False)


if __name__ == '__main__':
    """
    Example CLI usage
    $ serializeDir.py some_dir
    """
    input_dir = sys.argv[1]
    serializer = DirSerializer(input_dir)
    print(json.dumps(serializer.data, indent = 4))

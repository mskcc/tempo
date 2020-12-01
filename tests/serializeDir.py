#!/usr/bin/env python3
"""
Module for serializing directory file lists for use with testing
"""
import os
import sys
import hashlib
import json
import gzip
from pathlib import Path

class DirSerializer(object):
    """
    Serializes the listing of all files in a directory tree

    NOTE: some archive formats such as .zip, .gz, can have embedded timestamps that alter size and md5; also log files and html reports with timestamps can alter size and md5; might need to exclude these from checks for md5 and byte size
    https://unix.stackexchange.com/questions/31008/why-does-the-gzip-version-of-files-produce-a-different-md5-checksum
    """
    def __init__(self,
        root, # path to directoy to start from
        exclude_dirs = ['work', '.nextflow'], # directory names to skip over e.g. Nextflow work dir, etc
        exclude_exts = (), # ['.log'], # file types to exclude from output data
        exclude_sizes = (), # see note about needing to exclude size recording for some file types
        exclude_md5s = (),
        gz_txt_exts = ['.vcf.gz', '.txt.gz'] # recognize these text files in .gz format and count their lines
        ):
        self.root = root
        self.exclude_dirs = exclude_dirs
        self.exclude_exts = exclude_exts
        self.exclude_sizes = exclude_sizes
        self.exclude_md5s = exclude_md5s
        self.gz_txt_exts = gz_txt_exts

        self.data = {}

        if os.path.isdir(root):
            self.serialize_dir(root)
        else:
            self.serialize_file(root)

    def serialize_file(self, file):
        """
        Serialize a single file
        """
        path = os.path.abspath(file)
        dirname = os.path.dirname(file)
        key = file
        self.create_attributes(file, path, key)

    def serialize_dir(self, root):
        """
        Walk through the directory tree and serialize each file found
        """
        for dirname, dirs, files in os.walk(root):
            # remove excluded dirs and files from the list
            if self.exclude_dirs:
                dirs[:] = [ d for d in dirs if d not in self.exclude_dirs ]
            if self.exclude_exts:
                files[:] = [ fn for fn in files if any(not fn.endswith(ext) for ext in self.exclude_exts) ]

            # get attributes for each file
            for file in files:
                path = os.path.join(dirname, file)
                key = str(Path(path).relative_to(root)) # dont include the input dir name in the key
                self.create_attributes(file, path, key)

    def create_attributes(self, file, path, key):
        """
        Get all the attributes for a file; md5, byte size, line count
        Put the attributes in the instance `data` dict
        """
        self.data[key] = {}

        if self.exclude_md5s:
            if not any(file.endswith(ext) for ext in self.exclude_md5s):
                self.data[key]['md5'] = self.md5sum(path)
        else:
            self.data[key]['md5'] = self.md5sum(path)

        if self.exclude_sizes:
            if not any(file.endswith(ext) for ext in self.exclude_sizes):
                self.data[key]['size'] = os.path.getsize(path)
        else:
            self.data[key]['size'] = os.path.getsize(path)

        if self.isText(path):
            self.data[key]['lines'] = self.numLines(path)
        else:
            if any(file.endswith(ext) for ext in self.gz_txt_exts):
                self.data[key]['lines'] = self.numLinesGz(path)

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

    def numLinesGz(self, filepath):
        count = 0
        with gzip.open(filepath, 'rb') as f:
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
    serializer = DirSerializer(input_dir) # , exclude_md5s = ['.html', '.json']
    print(serializer.json(indent = 4))

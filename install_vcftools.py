#! /usr/bin/env python
import subprocess
import os

subprocess.call(["wget", "https://github.com/vcftools/vcftools/archive/master.zip"])
subprocess.call(["unzip", "master.zip"])
path = os.path.join(os.getcwd(), "vcftools-master")
perlpath = os.path.join(path, "src/perl/")
subprocess.call(["export", "PERL5LIB=" + perlpath], shell=True)
subprocess.check_call("./autogen.sh", shell=True, cwd=path)
subprocess.call("./configure", shell=True, cwd=path)
subprocess.call("make", shell=True, cwd=path)
subprocess.call("make install", shell=True, cwd=path) # needs sudo
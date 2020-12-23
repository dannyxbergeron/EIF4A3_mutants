# comes from http://thelazylog.com/install-python-as-local-user-on-linux/
# download and install python3.6.9 local in the cluster (the only one that worked
# with majiq) and then create a virtual env and install majic (see majiq.smk)

wget https://www.python.org/ftp/python/3.6.9/Python-3.6.9.tar.xz

tar -xvzf Python-3.6.9.tar.xz

find ./Python-3.6.9 -type d | xargs chmod 0755

cd Python-3.6.9

./configure --prefix=$HOME/python
make && make install



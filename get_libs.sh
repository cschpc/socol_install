mkdir -p downloads
pushd downloads 

curl -LO https://code.mpimet.mpg.de/attachments/download/28013/cdo-2.2.0.tar.gz && \
  tar xfvz cdo-2.2.0.tar.gz && \

curl -LO https://gitlab.dkrz.de/mpim-sw/libcdi/-/archive/cdi-1.7.2/libcdi-cdi-1.7.2.tar.gz && \
  tar xfvz libcdi-cdi-1.7.2.tar.gz  

curl -LO https://swprojects.dkrz.de/redmine/attachments/download/506/yaxt-0.9.1.tar.gz && \
  tar xfvz yaxt-0.9.1.tar.gz 

popd

# syntax=dockers/dockerfile:1

FROM python:3.7-popos
WORKDIR /code
RUN pip install collections
RUN pip install bs4
RUN pip install MdAnalysis==2.3.2
RUN pip install https://github.com/volkamerlab/opencadd/archive/master.tar.gz
RUN pip install rdkit==2.2
RUN pip install biotite
#RUN pip install -r rerquirements.txt
#COPY requirements.txt
ADD query-the-pdb.py
CMD ["python", "./query-the-pdb.py"]



function SF=space_frequency(X)
X=double(X);
[n0,n1]=size(X);%%%%   X��һ������
X=double(X);                          %�ռ�Ƶ��;
RF=0;
CF=0;

for fi=1:n0
    for fj=2:n1
        RF=RF+(X(fi,fj)-X(fi,fj-1)).^2;
    end
end

RF=RF/(n0*n1);

for fj=1:n1
    for fi=2:n0
        CF=CF+(X(fi,fj)-X(fi-1,fj)).^2;
    end
end

CF=CF/(n0*n0);%%%%����˼�����ռ�Ƶ���ǲ���ֻ����һ������

SF=sqrt(RF+CF)
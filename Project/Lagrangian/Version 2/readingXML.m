clear
clc
tree = xmlread('results.xml');
elemList = tree.getElementsByTagName('Step');
n = elemList.getLength;
for i=0:n-1
    x=elemList.item(i);
    value = x.getAttribute('type');
    if(strcmp(value,'dynamic'))
        javaString = x.getTextContent;
        matlabString = char(javaString);
        strData = strsplit(matlabString);
        data = str2double(strData);
        data = data(2:74);
        result(i+1).t = data(1);
        result(i+1).r1 = [data(2);data(3);data(4);];
        result(i+1).r2 = [data(8);data(9);data(10)];
        result(i+1).r3 = [data(14);data(15);data(16)];
    end
end
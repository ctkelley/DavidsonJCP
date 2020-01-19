% TDAV
%
% Run this script to generate the table in Section 3.1
%
gtab=[1000-1, 2000-1, 4000-1, 8000-1, 16000-1];
ihist=zeros(6,5);
for itab=1:5
     ihist(:,itab)=test_dav(gtab(itab));
end
headers={'h=1/1000', 'h=1/2000', 'h=1/4000', 'h=1/8000', 'h=1/16000'};
formats=' %9.4e & %9.4e & %9.4e & %9.4e & %9.4e';
fprint_tex(headers,formats,ihist); 
ihist
end



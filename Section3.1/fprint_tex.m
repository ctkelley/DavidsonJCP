function fprint_tex(headers,formats,data)
% FPRINT_TEX
% Print a LaTeX table from a matlab array.
%
% Inputs: 
% headers: the titles for the columns
%          example: headres={'foo', 'bar'}
% formats: c-style formatting for the columns. 
%          fprint_tex will add the carriage returns for you.
%          example: formats='%d & %7.2e';
%
% data: an mr x mc matlab array
%
[mr,mc]=size(data);
%
% Top matter goes here.
%
fprintf('\\begin{tabular}{');
for i=1:mc
    fprintf('l');
end
fprintf('} \n');
for i=1:mc-1
   fprintf('%9s &',headers{i});
end
fprintf('%9s \\\\ \n' ,headers{mc});
fprintf('\\hline \n');
%
% Put the guts of the table here.
%
for i=1:mr
    fprintf(strcat(formats,'\\\\ \n'), data(i,:))
end
%
% Print bottom matter and close out.
%
fprintf('\\hline \n');
fprintf('\\end{tabular} \n');




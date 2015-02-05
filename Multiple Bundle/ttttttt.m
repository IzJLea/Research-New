r=7;

for i=1:r
eval(['A' num2str(i) '= i']);
eval(['B' num2str(i) '= A' num2str(i) '*i']); % Issue was Here. Notice the positions of '
end
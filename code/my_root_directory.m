function y = my_root_directory

x = dir('/');
s = {x(:).name};
if any(strcmp(s, 'mindhive'))
    y = '/mindhive/nklab/u/svnh';
else
    y = '/Users/svnh2/Desktop/projects';
end
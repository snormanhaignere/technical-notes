root_directory = my_root_directory;
addpath([root_directory '/general-audio-code']);
addpath([root_directory '/general-analysis-code']);

x = [ones(5,1); zeros(200,1)];
x = x/sum(x);

n_steps = 5;
y = x;
for i = 1:n_steps
    y = myconv(x, y, 'causal', true);
    y = y/sum(y);
end
plot([x,y]);
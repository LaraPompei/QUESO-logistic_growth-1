clear all;
a_data
model_data

x = a_data;
y = model_data;

plot(x,y)
set(gca);
print -dpng base_model.png
pause;
clf;

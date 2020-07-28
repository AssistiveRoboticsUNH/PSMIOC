function plotTable(table)
x = 1:size(table,2);
data = [];%37.6 24.5 14.6 18.1 19.5 8.1 28.5 7.9 3.3 4.1 7.9 1.9 4.3]';
high = [];%2.1 4.4 0.4 3.3 2.5 0.4 1.6 0.8 0.6 0.8 2.2 0.9 1.5];
low  = [];%4.4 2.4 2.3 0.5 1.6 1.5 4.5 1.5 0.4 1.2 1.3 0.8 1.9];
j = 1;
for i = 1:size(table,2)
   j = mod(i-1,4)+1;
   data = [data;mean(table(j,i,:))];
   high = [high;max(table(j,i,:))];
   low = [low;min(table(j,i,:))];
end

bar(x,data)                

hold on

er = errorbar(x,data,low-data,high-data);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off
set(gcf,'Color','w')
end


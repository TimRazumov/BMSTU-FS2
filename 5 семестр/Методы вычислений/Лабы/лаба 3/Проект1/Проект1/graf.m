Pn = fopen ('polynom.txt','rt');
if Pn == -1 
    error('File is not opened'); 
end 

x=0;                % инициализаци€ переменной 
fpn=0;
f=0;
cnt=1;              % инициализаци€ счетчика 
while ~feof(Pn)    % цикл, пока не достигнут конец файла 
    [Vx,Nx] = fscanf(Pn, '%f',1);  %считывание одного 
    [Vf,Nf] = fscanf(Pn, '%f',1);
    [Vf_true,Nf_ftrue] = fscanf(Pn, '%f',1);
% значени€ double (V содержит значение 
% элемента, N Ц число считанных элементов) 
    if (Nf_ftrue > 0)        % если элемент был прочитан успешно, то 
            x(cnt)=Vx;   % формируем вектор-строку из значений V
            fpn(cnt)=Vf;
            f(cnt)=Vf_true;
            cnt=cnt+1;  % увеличиваем счетчик на 1     
    end 
end 
fclose(Pn); 

plot(x,fpn,'ro-', x, f,'b')% красное полином, синее исходна€ функци€
 legend('»нтерпол€нт','»сходна€ функци€');
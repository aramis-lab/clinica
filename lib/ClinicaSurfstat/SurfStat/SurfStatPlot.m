function [ t, df, pval ] = SurfStatPlot( x, y, M, g, varargin );

%Numeric or string variables against numeric, adjusted and/or by groups.
% 
% Usage: [ t, df, pval ] = SurfStatPlot( x, y [,M [,g [,varargin]]] );
% 
% x = n x 1 vector of numbers or cell array of strings or a term of these.
% y = n x 1 vector of numbers.
% M = term or anything that can be converted to a term, for adjusting y, 
%     default is M=1.
% g = n x 1 cell array of strings or a term of this. Separate plots are 
%     made for each different string (group) in g. Default is g=1.
% varargin = extra arguments for the plot, as for matlab's plot, e.g.
%     SurfStatPlot( x, y,  M,  g, 'LineWidth',2, 'MarkerSize',12)
% or  SurfStatPlot( x, y, [], [], 'LineWidth',2, 'MarkerSize',12)
% if M and g are not defined.
%
%   If g=1, y is adjusted for M in the model M+X, where M=1+term(M), 
% X=term(x), and tested against the null model M. 
% If g=/=1, y is adjusted for M in the model M+X+G+X*G, where 
% G=term(g), and tested against the null model M+X+G.  
% t    = T statistic if x is one variable and g=1, otherwise F statistic. 
% df   = degrees of freedom.
% pval = P-value. 

y=double(y);

if nargin<3 | isempty(M)
    M=1;
else
    if isnumeric(M)
        M=1+term(M);
    end
end
if nargin<4 | isempty(g)
    G=1;
else
    G=term(g);
end

X=term(x);
XG=X+G+X*G;
slm=SurfStatLinMod(y,M+XG);
Yhat=slm.X*slm.coef;
XGmat=double(XG);
Yhatadj=XGmat*pinv(XGmat)*Yhat;
c=mean(y)-mean(Yhatadj);
Yhatadj=Yhatadj+c;
Yadj=Yhatadj+y-Yhat;

nx=size(X,2);
ng=size(G,2);
x=double(X)*((1:nx)');
g=double(G)*((1:ng)');

if nx==1
    if ng==1
        iu=[find(x==min(x)); find(x==max(x))];
        plot(x(iu),Yhatadj(iu),'-r',x,Yadj,'.b',varargin{:});
        slm=SurfStatT(slm,x);
        if isfield(slm,'dfs')
            slm.df=slm.dfs;
        end
        pval=stat_threshold(0,1,0,slm.df,[10 abs(slm.t)],[],[],[],[],[],[],0);
        pval=pval(2)*2;
        xlabel([inputname(1) ...
            ': T = ' num2str(round(slm.t*100)/100) ...
            ', df = ' num2str(round(slm.df*10)/10) ...
            ', P = ' num2str(round(pval*1000)/1000)]);
    else 
        hs=zeros(1,ng);
        for j=1:ng
            i=(g==j);
            xi=x(i);
            Yhatadji=Yhatadj(i);
            iu=[find(xi==min(xi)); find(xi==max(xi))];
            h=plot(xi(iu),Yhatadji(iu),'-',xi,Yadj(i),'.',varargin{:});
            set(h(2),'Color',get(h(1),'Color'));
            hs(j)=h(1);
            if j==1
                hold all;
            end
        end
        hold off;
        legend(hs,char(G));
        slm0=SurfStatLinMod(y,M+X+G);
        slm=SurfStatF(slm,slm0);
        pval=stat_threshold(0,1,0,slm.df,[10 slm.t],[],[],[],[],[],[],0);
        pval=pval(2);
        xlabel([inputname(1) ' * ' inputname(4) ...
            ': F = ' num2str(round(slm.t*100)/100) ...
            ', df = ' num2str(slm.df(1)) ',' num2str(slm.df(2)) ...
            ', P = ' num2str(round(pval*1000)/1000)]);
    end 
else
    if ng==1
        [xu,iu]=unique(x);
        plot(xu,Yhatadj(iu),'-r',x,Yadj,'.b',varargin{:});
        slm0=SurfStatLinMod(y,M);
        str='';
    else
        hs=zeros(1,ng);
        for j=1:ng
            i=(g==j);
            xi=x(i);
            Yhatadji=Yhatadj(i);
            [xu,iu]=unique(xi);
            h=plot(xu,Yhatadji(iu),'-',xi,Yadj(i),'.',varargin{:});
            set(h(2),'Color',get(h(1),'Color'));
            hs(j)=h(1);
            if j==1
                hold all;
            end
        end
        hold off;
        legend(hs,char(G));
        slm0=SurfStatLinMod(y,M+X+G);
        str=[' * ' inputname(4)];
    end
    xlim([0.5 nx+0.5]);
    set(gca,'XTick',(1:nx)');
    set(gca,'XTickLabel',strvcat(char(X)));
    slm=SurfStatF(slm,slm0);
    pval=stat_threshold(0,1,0,slm.df,[10 slm.t],[],[],[],[],[],[],0);
    pval=pval(2);
    xlabel([inputname(1) str ...
        ': F = ' num2str(round(slm.t*100)/100) ...
      ', df = ' num2str(slm.df(1)) ',' num2str(slm.df(2)) ...
      ', P = ' num2str(round(pval*1000)/1000)]);    
end    
if nargin<3 | prod(size(M))==1
    ylabel(inputname(2));
else
    ylabel([inputname(2) ' adjusted for ' inputname(3)]);
end
t=slm.t;
df=slm.df;
set(gcf,'PaperPosition',[0.25 2.5 6 4.5]);

return
end

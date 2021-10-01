% Generalized Generalized Mittag-Leffler function
% Reference:
% [1] A.K. Shukla ?, J.C. Prajapati. "On a generalization of Mittag-Leffler 
% function and its properties" J. Math. Anal. Appl. 336 (2007) 797–811
% [2] A. A. KILBAS, M. SAIGOb, and R. K. SAXENA, “Generalized 
% Mittag-Leffler function and generalized fractional calculus operators,” 
% Integral Transforms and Special Functions, v. 15, No. 1, 2004, pp. 31–49.
% Prepared by YangQuan Chen. 09/16/2008
% Based on the ml_fun code in Dingyu Xue, YangQuan Chen* and Derek Atherton.
% “Linear Feedback Control – Analysis and Design with Matlab”. SIAM Press, 2007, 
% ISBN: 978-0-898716-38-2. (348 pages)
% Using notations in [1]
% a: \alpha \in C, and \re{\alpha) >0
% b: \beta \in C, and \re{\beta) >0
% c: \gamma \in C, and \re{\gamma) >0
% d: q \in (0,1) \union N
% x: the variable
% eps0: specified accuracy
% Definition in [1]
% E_{\alpha,\beta}^{\gamma, q} (z) =
% \sum_{n=0}^{\infty} (z^n/n!) (\gamma)_{qn} / \Gamma(\alpha n +\beta)
% see also:
% http://www.mathworks.com/
% matlabcentral/fileexchange/loadFile.do?objectId=20849&objectType=FILE
% GML_FUN is the special case of GGML_FUN when q=1
% Note: GGML_FUN converges absolutely when q< \re{\alpha) + 1
%
 function f=ggml_fun(a,b,c,d,x,eps0)
 gamma_c=1.0/gamma(c);
 if nargin<6, eps0=eps; end
 if (d>= real(a)+1), sprintf('%s','Note: GGML_FUN converges absolutely when q< \re{\alpha) + 1 '),   f=NaN,  return, end
 if ( real(a)<0), sprintf('%s',' \re{\alpha) must be greater than 0 '),  f=NaN,    return, end
 if ( real(b)<0), sprintf('%s',' \re{\beta) must be greater than 0 '),   f=NaN,   return, end
 if ( real(c)<0), sprintf('%s',' \re{\gamma) must be greater than 0 '),  f=NaN,    return, end
 f=0; fa=1; j=0;
 while norm(fa,1)>=eps0
    fa=(gamma(c+d*j)*gamma_c)/gamma(j+1)/gamma(a*j+b) *x.^j;
    f=f+fa; j=j+1; 
 end
% 
% GGML test code
% exp(x)
% deltat=0.01;
% x=[-1:deltat:1]; 
% y1=ggml_fun(1,1,1,1,x);
% figure;plot(x,y1,'r',x,exp(x),'k')
% figure;plot(x,y1-exp(x),'k')
% % ok, now try MLF when c=1
% y2=ggml_fun(0.5,0.5,1,1,x);
% y3=MLF_M(0.5,0.5,x); % Igor Podlubny's code
% figure;plot(x,y1,'b',x,y2,'r',x,y3,'k')
% %figure;plot(x, y2-y3','k')
% legend('e^x=E_{1,1}^{(1)}', 'E_{0.5,0.5}^{(1)} - this code', 'E_{0.5,0.5}^{(1)} - Podlubny code')
% % okay, let us now try GGML
% y2=ggml_fun(1.5,0.5,1,1,x);
% y5=ggml_fun(1.5,0.5,0.5,1,x);
% y6=ggml_fun(1.5,0.5,1.5,1,x);
% y2q=ggml_fun(1.5,0.5,1,.5,x);
% y5q=ggml_fun(1.5,0.5,0.5,.5,x);
% y6q=ggml_fun(1.5,0.5,1.5,.5,x);
% y2q2=ggml_fun(1.5,0.5,1,2,x);
% y5q2=ggml_fun(1.5,0.5,0.5,2,x);
% y6q2=ggml_fun(1.5,0.5,1.5,2,x);
% 
% figure;plot(x,y2,'r',x,y2q,'r:',x,y2q2,'r-.',    x,y5,'k',x,y5q,'k:',x,y5q2,'k-.',    x,y6,'b',x,y6q,'b:',x,y6q2,'b-.')
% legend('GGML \gamma=1, q=1','GGML \gamma=1, q=0.5','GGML \gamma=1, q=2', ...
%     'GGML \gamma=.5 q=1','GGML \gamma=.5 q=0.5','GGML \gamma=.5 q=2',...
%     'GGML \gamma=1.5, q=1','GGML \gamma=1.5, q=0.5','GGML \gamma=1.5, q=2')
% title('E_{1.5, 0.5}^{\gamma, q}(x)')


% mapper 和 demapper 在最下面

% p3的模擬
n=[0,1,2,3,4,5,6,7,8,9,10];
M=4;
%eb=(M^2-1)/(12*log2(M));%PAM
%eb=(M-1)/(6*log2(M));%QAM
eb=(1/4)/(log2(M)*sin(pi/M)*sin(pi/M));%PSK
N0=[];
for i=0:10
    N0=[N0,[eb*10^(-i/10)]];
end
SER=[];
for i=1:11
    SER=[SER,[run_PSK(N0(i),M,1000000)]];
    
end
theoretical=[];
for i=1:11
    %x=2*Q(sqrt((6*log(M)*eb)/(log(2)*(M^2-1)*N0(i))));%PAM
    %x=4*Q(sqrt((3*log(M)*eb)/(log(2)*(M-1)*N0(i))));%QAM
    x=2*Q(sqrt(2*log2(M)*sin(pi/M)*sin(pi/M)*eb/N0(i)));%PSK
    theoretical=[theoretical,[x]];
end
figure(1)
semilogy(n,theoretical,'-','Color',[0 0 1])
hold on
semilogy(n,SER,'o','Color',[0.8 0.6 0])
hold off
%P4寫不完...
function res=Q(x)
   phi=(1/2)*(1-erf(-x/sqrt(2)));
   res=(1-phi);
end


function sym=noise(sym,N)
   for i=1:length(sym)
       n=sqrt(N/2)*(randn(1,1)+j*randn(1,1));
       sym(i)=sym(i)+n;
   end
end
function sym_seq=gray_code(l)
   sym_seq=["0","1"];
   for i=1:l-1
       a=sym_seq;
       b=flip(sym_seq);
       for j=1:length(a)
           a(j)="0 "+a(j);
           b(j)="1 "+b(j);
       sym_seq=[a,b];
       end
       
   end
   
end
function res=PAM(bin_seq,M,d)
   res=[];
   ep=d/2;
   seg_len=log(M)/log(2);
   g_c=gray_code(seg_len);
   for i=1:seg_len:length(bin_seq)-seg_len+1
       %check=i
       b=bin_seq(i:i+seg_len-1);
       for j=1:length(g_c)
           symbol=(split(g_c(j)," "))';
           if b==symbol
               res=[res,[(2*j-1-M)*ep]];
           end
           
       end
       
       
   end
   
end
function res=PSK(bin_seq,M,d)
   ep=d/(2*sin(pi/M));
   res=[];
   seg_len=log(M)/log(2);
   g_c=gray_code(seg_len);
   for i=1:seg_len:length(bin_seq)-seg_len+1
      b=bin_seq(i:i+seg_len-1);
      for k=1:length(g_c)
          symbol=(split(g_c(k)," "))';
          if b==symbol
              res=[res,[(cos(2*pi*(k-1)/M)+j*sin(2*pi*(k-1)/M))*ep]];
          end
      end
   end
end
function res=QAM(bin_seq,M,d)
   %rest=gray_code(4)
   Am=[];
   len=sqrt(M);
   count=1;
   ep=d/(2);
   for x=1:len
       if count==1
           for y=1:len
               Am=[Am,[(2*x-1-len)*ep+j*(2*y-1-len)*ep]];
           end
           count=0;
       else
           for y=len:-1:1
               Am=[Am,[(2*x-1-len)*ep+j*(2*y-1-len)*ep]];
           end
           count=1;
           
       end
   end
   
   
   res=[];
   seg_len=log(M)/log(2);
   g_c=gray_code(seg_len);
   for i=1:seg_len:length(bin_seq)-seg_len+1
      b=bin_seq(i:i+seg_len-1);
      for k=1:length(g_c)
          symbol=(split(g_c(k)," "))';
          if b==symbol
              res=[res,[Am(k)]];
          end
      end
   end
   
end
function res=dPAM(sym_seq,M,d)
   res=[];
   g_c=gray_code(log(M)/log(2));
   rn=[];
   for k=1:M
       rn=[rn,[(2*k-1-M)*(d/2)]];
   end
   for k=1:length(sym_seq)
       min=100000000;
       n=-1;
       for i=1:length(rn)
           distance=abs(sym_seq(k)-rn(i));
           if distance<min
               n=i;
               min=distance  ;   
           end
       end
       res=[res,split(g_c(n)," ")'];
   end
   
end
function res=dPSK(sym_seq,M,d)
   res=[];
   g_c=gray_code(log(M)/log(2));
   rn=[];
   for k=1:M
       rn=[rn,[(d/2)*cos(2*pi*(k-1)/M)+j*(d/2)*sin(2*pi*(k-1)/M)]];
   end
   for k=1:length(sym_seq)
       min=100000000;
       n=-1;
       for i=1:length(rn)
           distance=abs(sym_seq(k)-rn(i));
           if distance<min
               n=i;
               min=distance  ;   
           end
       end
       res=[res,split(g_c(n)," ")'];
   end
   
end
function res=dQAM(sym_seq,M,d)
   res=[];
   g_c=gray_code(log(M)/log(2));
   Am=[];
   len=sqrt(M);
   count=1;
   ep=d/(2);
   for x=1:len
       if count==1
           for y=1:len
               Am=[Am,[(2*x-1-len)*ep+j*(2*y-1-len)*ep]];
           end
           count=0;
       else
           for y=len:-1:1
               Am=[Am,[(2*x-1-len)*ep+j*(2*y-1-len)*ep]];
           end
           count=1;
           
       end
    end
  
   for k=1:length(sym_seq)
       min=100000000;
       n=-1;
       for i=1:length(Am)
           distance=abs(sym_seq(k)-Am(i));
           if distance<min
               n=i;
               min=distance ;    
           end
       end
       res=[res,split(g_c(n)," ")'];
   end
   
end
function res=random_string(n)
   res=[];
   for i=1:n
       x=rand;
       if x>1/2
           res=[res,["0"]];
       else
           res=[res,["1"]];
       end
   end
end
function P=run_PAM(N,M,n)
   res=0;
   bin_r=random_string(n);
   bin=PAM(bin_r,M,1);
   bin=noise(bin,N);
   dmap=dPAM(bin,M,1);
   len=log(M)/log(2);
   for i=1:len:length(dmap)-len+1
       if bin_r(i:i+len-1)~=dmap(i:i+len-1)
           res=res+1;
       end
   end
   P=res/n;
end
function P=run_QAM(N,M,n)
   res=0;
   bin_r=random_string(n);
   bin=QAM(bin_r,M,1);
   bin=noise(bin,N);
   dmap=dQAM(bin,M,1);
   len=log(M)/log(2);
   for i=1:len:length(dmap)-len+1
       if bin_r(i:i+len-1)~=dmap(i:i+len-1)
           res=res+1;
       end
   end
   P=res/n;
end
function P=run_PSK(N,M,n)
   res=0;
   bin_r=random_string(n);
   bin=PSK(bin_r,M,1);
   bin=noise(bin,N); 
   dmap=dPSK(bin,M,1);
   len=log(M)/log(2);
   for i=1:len:length(dmap)-len+1
       if bin_r(i:i+len-1)~=dmap(i:i+len-1)
           res=res+1;
       end
   end
   P=res/n;
end
function res=symbol_mapper(bin_seq,M,d,name)
   if name=="PAM"
       res=PAM(bin_seq,M,d)
   elseif name=="PSK"
       res=PSK(bin_seq,M,d)
   else
       res=QAM(bin_seq,M.d)
   end
          
end
function res=MD_symbol_demapper(sym_seq,M,d,name)
   if name=="PAM"
       res=dPAM(sym_seq,M,d)
   elseif name=="PSK"
       res=dPSK(sym_seq,M,d)
   else
       res=dQAM(sym_seq,M.d)
   end
          
end
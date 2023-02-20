s=["a","b","c","d","e","f","g","h"]
p=[0.2,0.05,0.005,0.2,0.3,0.05,0.045,0.15]
r=[10,20,50,100,200,500,1000]
p1_e=["g","a","c","a","b"]
h=0
dict=huffman_dict(s,p)
L_R=[]
len=0
for i=1:20
    
    rs=randomstring(s,p,10)
    c=huffman_enc(rs,dict)
    len=len+length(split(c))-1
    L_R=[L_R,[len/i]]
end
check=len
histogram(L_R)
for i=1:length(p)
    h=h-p(i)*log(p(i))/log(exp(1))
end
n_10=[]
n_h=[]
n_h2=[]

for i=1:length(r)
    l=repeat(s,p,10,r(i),dict)
    l=l/10
    n_10=[n_10,[l]]
    n_h=[n_h,[h]]
    n_h2=[n_h2,[h+1]]
end
n_50=[]
for i=1:length(r)
    l=repeat(s,p,50,r(i),dict)
    l=l/50
    n_50=[n_50,[l]]
end
n_100=[]
for i=1:length(r)
    l=repeat(s,p,100,r(i),dict)
    l=l/100
    n_100=[n_100,[l]]
end
semilogx(r,n_10,r,n_50,r,n_100,r,n_h,r,n_h2)
dict=huffman_dict(s,p)
 bin=huffman_enc(p1_e,dict)
 sys=huffman_dec(bin,dict)
 %rs=randomstring(s,p,10)
% a_3=huffman_enc(rs,dict)
% b_3=repeat(s,p,10,1000)
% c_3=(b_3)/10
function dict=huffman_dict(symbols,prob)
    dict={}
    
    for v=1:length(prob)
        dict{v,1}=symbols(v)
        dict{v,2}=prob(v)
        dict{v,5}=""
    end
    x=0
    y=0
    i=1
    j=2
    rec=["error"]
    while x+y<1
        p=sort(prob)
        min=p(i)
        min2=p(j)
        x=min
        y=min2
        i=i+2
        j=j+2
        p3=[min+min2]
        prob=[prob,p3]
        %l1=dict(:,1)
        %a=find(l1==min)
        count1=0
        count2=0
        k=1
        while count1+count2<2
            if prob(k)==min & count1==0 & ~ismember(symbols(k),rec)
                a=symbols(k)
                i1=k
                r=[a]
                rec=[rec,r]
                count1=count1+1
            elseif prob(k)==min2 & count2==0 & ~ismember(symbols(k),rec)
                b=symbols(k)
                i2=k
                r=[b]
                rec=[rec,r]
                count2=count2+1
            end
            k=k+1
        end
        
        ele=length(dict(:,1))+1
        dict{ele,1}=a+b
        dict{ele,2}=min+min2
        dict{ele,3}=i1
        dict{ele,4}=i2
        dict{ele,5}=""
        s3=[a+b]
        symbols=[symbols,s3]
        s=dict(:,1)
        len=length(s)
        for v=len:-1:1
            l=dict{v,3}
            r=dict{v,4}
            if ~isempty(l) & ~isempty(r)
               dict{l,5}=dict{v,5}+"0"+" "
               dict{r,5}=dict{v,5}+"1"+" "
            end
        end
     end
    
end      

function bin=huffman_enc(sys,dict)
   bin=""
   s=[]
   for i=1:length(dict(:,1))
       s=[s,[dict{i,1}]]
   end
   for i=1:length(sys)
       for j=1:length(s)
           if sys(i)==s(j)
               bin=bin+dict{j,5}
           end
       end
       
   end
end
function sys2=huffman_dec(bin,dict)
   code=dict(:,1)
   c=[]
   for i=1:length(code)
       if isempty(dict{i,3}) & isempty(dict{i,4})
           c=[c,[dict{i,5}]]
       end
   end
   
   bin=bin.split()
   initial=1
   s=""
   sys=[]
   for i=1:length(bin)
       s=s+bin(i)+" "
       for j=1:length(c)
           if s==c(j)
               a=c(j)
               sys=[sys,[a]]
               s=""
           end
       end
       
   end
   sys2=[]
   for i=1:length(sys)
       for j=1:length(c)
           if c(j)==sys(i)
               sys2=[sys2,[dict{j,1}]]
           end
       end
   end
  
   
end
function rs=randomstring(s,p,n)
  for i=2:length(p)
      p(i)=p(i-1)+p(i)
  end
  p=[[0],p]
  rs=[]
  for i=1:n
      x=rand
      for j=2:length(p)
          if x<=p(j) & x>=p(j-1)
              rs=[rs,[s(j-1)]]
          end
      end
  end
end
function len=repeat(s,p,n,t,dict)
   %dict=huffman_dict(s,p)
   len=0
   for i=1:t
       rs=randomstring(s,p,n)
       c=huffman_enc(rs,dict)
       len=len+length(split(c))-1
       
   end
   len=len/t
end

    
    
          
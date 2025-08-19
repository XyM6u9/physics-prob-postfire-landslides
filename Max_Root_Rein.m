function[Pr_Max]=Max_Root_Rein(kappaR,lambdaR,Wo,meanF_u,meanF_w,meanF_l,Fa,Fb,mF,kF,muF,sigF)
string choice;
isittrue=0;
while isittrue==0
    choice=inputdlg("Enter the sahpe of distribution function for root diameters Choices are uniform, weibull, or lognormal:",'s');
    if choice=="uniform" || choice=="weibull" || choice=="lognormal"
        isittrue=1;
    end
end

if choice=="uniform"
    if kappaR==1
        Pr_Max=meanF_u;
    elseif kappaR==-1/lambdaR
        Pr_Max=exp(((1+lambdaR)/lmbda*log(Fb)-1))/(1+lambdaR)/(Fb^(1/lambdaR)-Fa^(1/lambdaR));
    elseif kappaR<1
        Pr_Max=(Fb^((1+lambdaR)/lambdaR))/(Fb^(1/lambdaR)-Fa^(1/lambdaR))*(lambdaR*(1-kappaR))^((lambdaR*(1-kappaR))/(1+lambdaR*kappaR))/...
            ((1+lambdaR)^((1+lambdaR)/(1+lambdaR*kappaR)));
    else
        Pr_Max=(Fb^(1-kappaR))/(1+lambdaR*kappaR)*(Fb^((1+lambdaR*kappaR)/lambdaR)-Fa^((1+lambdaR*kappaR)/lambdaR))/(Fb^(1/lambdaR)-Fa^(1/lambdaR));
    end
elseif choice=="weibull"
    if kappaR==1
        Pr_Max=meanF_w;
    elseif kappaR<1
        i=1;
        for x=0:.001:0.50
            Pr_disp(i,1)=Wo*(kF^kappaR)*x*gammainc(((Wo*x)^(mF/(1-kappaR)))/(kF^mF),(kappaR/mF)+1,'upper');
            Pr_disp(i,2)=x;
            i=i+1;
            Pr_Max=max(Pr_disp(:,1));
        end
    else
        i=1;
        for x=0:0.001:0.50
            Pr_disp(i,1)=Wo*(kF^kappaR)*x*gammainc(((Wo*x)^(mF/(1-kappaR)))/(kF^mF),(kappaR/mF)+1,'lower');
            Pr_disp(i,2)=x;
            i=i+1;
            Pr_Max=max(Pr_disp(:,1));            
        end   
    end
else
    if kappaR==1
        Pr_Max=meanF_l;
    elseif kappaR<1
        i=1;
        for x=0:0.001:0.50
          Pr_disp(i,1)=0.5*Wo*x*exp(0.5*(kappaR^2)*(sigF^2)+kappaR*muF)*(1+erf((kappaR*(sigF^2)+muF-log(((Wo*x))^(1/(1-kappaR))))/(sigF*sqrt(2)))); 
          Pr_disp(i,2)=x;
          i=i+1;
          Pr_Max=max(Pr_disp(:,1));
        end
    else
        i=1;
        for x=0:0.001:0.50
          Pr_disp(i,1)=0.5*Wo*x*exp(0.5*(kappaR^2)*(sigF^2)+kappaR*muF)*(1-erf((kappaR*(sigF^2)+muF-log(((Wo*x))^(1/(1-kappaR))))/(sigF*sqrt(2)))); 
          Pr_disp(i,2)=x;
          i=i+1;
          Pr_Max=max(Pr_disp(:,1));
        end 
    end
end
end

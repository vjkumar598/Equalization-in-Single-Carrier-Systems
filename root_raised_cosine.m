%root_raise_cosine_ (RRC)filter in time domain.
%m : oversampling factor
%a : excess bandwidth
%l : truncate the pulse from -l to l

function [p,time] = root_raised_cosine(a,m,l)
                   l_os = floor(m * l); %length of one sided time vector : i.e #of samples taken on each side of peak
                   t_os = (1:l_os)./m; % one sided time vector
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %%
                   f1 = 4 * a * cos(pi * (1 + a) * t_os);
                   f2 = (sin(pi * (1-a) * t_os))./(t_os);
                   f3 = pi * (1 - (4 * a * t_os).^2);
                   %%
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   % zero test : 0/0 form at three points, i.e at
                   % m/(4*a),-m/(4*a) and 0
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   k = m /(4 * a); 
                   if(k == floor(k))
                       f1(k)= 4*a*pi*(1 + a)*sin(pi*(1+a)/(4*a)); %by L'opitals rule
                       f2(k)= 16*a^2*sin(pi*(1-a)/(4*a)) - 4*a*pi*(1-a)*cos(pi*(1-a)/(4*a));
                       f3(k)= 8*pi*a;
                   end
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   p0 = ((4*a)/pi)+(1-a) ;% pulse value at time 0
                                          %by L'opitals rule  
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %%
                   p_os = (f1 + f2)./ f3; %one sided pulse
                   p = [fliplr(p_os),p0,p_os];%two sided pulse
                   time = [fliplr(-t_os),0,t_os];  % alternative is  time = (-l_os:l_os)./m;
                    
end
                   
                   
               
           
           
           
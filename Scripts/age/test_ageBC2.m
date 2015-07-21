val=mask3.*0;
for i=2:402-1
for j=2:482-1
  if (mask3(i,j)>0) 
    % Add tracer if inflow through western cell wall
    if (mask3(i-1,j)-hmask(i-1,j)<0.0) 
      % First, check if neighboring borders are masked
      % West & North & South
      if ((mask3(i,j+1)-hmask(i,j+1))<0.0 & (mask3(i,j-1)-hmask(i,j-1))<0.0) 
        if (u(i-1,j  )> 0.0 | v(i  ,j  )< 0.0 | v(i  ,j-1)> 0.0) 
          ffac=1.0;
        else
          ffac=0.0;
        end
      % West & South & East
      elseif ((mask3(i,j-1)-hmask(i,j-1))<0.0 & (mask3(i+1,j)-hmask(i+1,j))<0.0) 
        if (u(i-1,j  )> 0.0 | u(i  ,j  )< 0.0 | v(i  ,j-1)> 0.0) 
          ffac=1.0;
        else
          ffac=0.0;
        end
      % West & North & East
      elseif ((mask3(i,j+1)-hmask(i,j+1))<0.0 & (mask3(i,j+1)-hmask(i,j+1))<0.0) 
        if (u(i-1,j  )> 0.0 | v(i  ,j  )< 0.0 | u(i  ,j  )< 0.0) 
          ffac=1.0;
        else
          ffac=0.0;
        end
      % West & North
      elseif (mask3(i,j+1)-hmask(i,j+1)<0.0) 
        if (u(i-1,j  )> 0.0 | v(i  ,j  )< 0.0) 
          ffac=1.0;
        else
          ffac=0.0;
        end
      % West & South
      elseif (mask3(i,j-1)-hmask(i,j-1)<0.0) 
        if (u(i-1,j  )> 0.0 | v(i  ,j-1)> 0.0) 
          ffac=1.0;
        else
          ffac=0.0;
        end	
      %  West & East  
      elseif (mask3(i+1,j)-hmask(i+1,j)<0.0) 
        if (u(i-1,j  )> 0.0 | u(i+1,j  )< 0.0) 
          ffac=1.0;;
        else
          ffac=0.0;
        end	
      % West only  
      else
        if (u(i-1,j  )> 0.0) 
          ffac=1.0;
        else
          ffac=0.0;
        end	  
      end
    %
    % Add tracer if inflow through eastern cell wall
    %
    elseif (mask3(i+1,j)-hmask(i+1,j)<0.0) 
      % East & North & South
      if ((mask3(i,j+1)-hmask(i,j+1))<0.0 & (mask3(i,j-1)-hmask(i,j-1))<0.0) 
        if (u(i  ,j  )< 0.0 | v(i  ,j  )< 0.0 | v(i  ,j-1)> 0.0) 
          ffac=1.0;
        else
          ffac=0.0;
        end
      % East & North
      elseif (mask3(i,j+1)-hmask(i,j+1)<0.0) 
        if (u(i  ,j)< 0.0 | v(i  ,j)< 0.0) 
          ffac=1.0;
        else
          ffac=0.0;
        end
      % East & South
      elseif (mask3(i,j-1)-hmask(i,j-1)<0.0) 
        if (u(i  ,j  )< 0.0 | v(i  ,j-1)> 0.0) 
          ffac=1.0;
        else
          ffac=0.0;
        end	
      % East only  
      else
        if (u(i  ,j  )< 0.0) 
          ffac=1.0;
        else
          ffac=0.0;
        end	  
      end
    %
    % Add tracer if inflow through northern cell wall
    %
    elseif (mask3(i,j+1)-hmask(i,j+1)<0.0) 
      % North & South
      if (mask3(i,j-1)-hmask(i,j-1)<0.0) 
        if (v(i  ,j  )< 0.0 | v(i  ,j-1)> 0.0) 
          ffac=1.0;
        else
          ffac=0.0;
        end	
      % North only  
      else
        if (v(i,j)< 0.0) 
          ffac=1.0;
        else
          ffac=0.0;
        end	  
      end
    %
    % Add tracer if inflow through southern cell wall
    %
    elseif (mask3(i,j-1)-hmask(i,j-1)<0.0) 
      if (v(i,j-1)> 0.0) 
        ffac=1.0;
      else
        ffac=0.0;
      end
    else
      ffac=0.0;
    end
  else
    ffac=0.0;
  end
  val(i,j)=ffac;
end
end

function mean_t = mean_norm_truncated(mean,stdv,xmin,xmax)

mean_t = 1/(ferum_cdf(1,(xmax-mean)/stdv,0,1)-ferum_cdf(1,(xmin-mean)/stdv,0,1)) * ...
            ( mean * (ferum_cdf(1,(xmax-mean)/stdv,0,1)-ferum_cdf(1,(xmin-mean)/stdv,0,1)) - ...
              stdv * (ferum_pdf(1,(xmax-mean)/stdv,0,1)-ferum_pdf(1,(xmin-mean)/stdv,0,1)) );
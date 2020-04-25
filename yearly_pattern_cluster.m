clear all

%util='e';
util='w';

load('building_info.mat')
w_data=load_data(util);
builds=fields(w_data);
e_data=load_data('e');
q_avg=zeros(length(builds),7);

for i=1:length(builds)
    
    b_loop=char(builds(i));
    
d=w_data.(b_loop).d;
q=w_data.(b_loop).q;

%the building ETC is missing for energy data so just make the whole
%timeseries nans
if strcmp(b_loop,'ETC')|strcmp(b_loop,'MAI')|strcmp(b_loop,'RHD')
    d_e=d;
    e=NaN(length(d_e),1);
else

d_e=e_data.(b_loop).d;
e=e_data.(b_loop).q;
end

%Cut to specific year
q=q(year(d)==2015);
d=d(year(d)==2015);

e=e(year(d_e)==2015);
d_e=d_e(year(d_e)==2015);


%  %by intensity
  
  q_i=q./building_info{i,5};
  e_i=e./building_info{i,5};

  q=normalize_data(q);
e=normalize_data(e);
  
  %Find the median intensity of each building to set parameters for
  %classification
  e_i_med(i)=median(e_i,'omitnan');
  q_i_med(i)=median(q_i,'omitnan');
  
[week_data]=timeseries_to_weekdata(d,q);
[week_data_e]=timeseries_to_weekdata(d_e,e);
[week_data_qi]=timeseries_to_weekdata(d,q_i);
[week_data_ei]=timeseries_to_weekdata(d_e,e_i);

%% Yearly Analysis Plot

plot_inds=find_pattern(week_data,'n');
plot_inds_e=find_pattern(week_data_e,'n');
plot_inds_qi=find_pattern(week_data_qi,'wi');
plot_inds_ei=find_pattern(week_data_ei,'ei');

%Keep track of how many overall each pattern occurred
plot_inds_all(i,:)=plot_inds;
plot_inds_e_all(i,:)=plot_inds_e;
plot_inds_qi_all(i,:)=plot_inds_qi;
plot_inds_ei_all(i,:)=plot_inds_ei;

%Create plot of histogram of each metric by week to see which ones were the
%most common across all of them


%%


%find the mode of the weekly pattern

%Normalized water
mode_week=mode(plot_inds);
week_data_med=week_data(plot_inds==mode_week,:);

%Normalized Energy
mode_week_e=mode(plot_inds_e);
week_data_med_e=week_data_e(plot_inds_e==mode_week_e,:);

%Water Intensity
mode_week_qi=mode(plot_inds_qi);
week_data_med_qi=week_data_qi(plot_inds_qi==mode_week_qi,:);

%Energy Intensity
mode_week_ei=mode(plot_inds_ei);
week_data_med_ei=week_data_ei(plot_inds_ei==mode_week_ei,:);


%find the median and save it in a matrix
q_avg(i,:)=median(week_data_med,'omitnan');
e_avg(i,:)=median(week_data_med_e,'omitnan');
q_avg_qi(i,:)=median(week_data_med_qi,'omitnan');
q_avg_ei(i,:)=median(week_data_med_ei,'omitnan');

end

%%
figure

subplot(2,2,1)
h=histogram(categorical(plot_inds_all))
vals=h.Values;
c = categorical({'Low,Flat','Low,Curve','Med,Flat','Med,Curve','High,Flat','High,Curve'});
c = reordercats(c,{'Low,Flat','Low,Curve','Med,Flat','Med,Curve','High,Flat','High,Curve'});
%prices = [1.23 0.99 2.3];
b=bar(c,vals)
ylabel('Number of Weeks')
title('Overall Normalized Water by Week')

subplot(2,2,3)
h=histogram(categorical(plot_inds_e_all))
vals=h.Values;
c = categorical({'Low,Flat','Low,Curve','Med,Flat','Med,Curve','High,Flat','High,Curve'});
c = reordercats(c,{'Low,Flat','Low,Curve','Med,Flat','Med,Curve','High,Flat','High,Curve'});
%prices = [1.23 0.99 2.3];
b=bar(c,vals)
ylabel('Number of Weeks')
title('Overall Normalized Energy by Week')

subplot(2,2,2)
h=histogram(categorical(plot_inds_qi_all))
vals=h.Values;
c = categorical({'Low,Flat','Low,Curve','Med,Flat','Med,Curve','High,Flat','High,Curve'});
c = reordercats(c,{'Low,Flat','Low,Curve','Med,Flat','Med,Curve','High,Flat','High,Curve'});
%prices = [1.23 0.99 2.3];
b=bar(c,vals)
ylabel('Number of Weeks')
title('Overall Water Intensity by Week')

subplot(2,2,4)
h=histogram(categorical(plot_inds_ei_all))
vals=h.Values;
vals(3:6)=vals(2:5);
vals(2)=0;
c = categorical({'Low,Flat','Low,Curve','Med,Flat','Med,Curve','High,Flat','High,Curve'});
c = reordercats(c,{'Low,Flat','Low,Curve','Med,Flat','Med,Curve','High,Flat','High,Curve'});
%prices = [1.23 0.99 2.3];
b=bar(c,vals)
ylabel('Number of Weeks')
title('Overall Energy Intensity by Week')

%%
%  y = quantile(e_i_med,[1/3, 2/3])
%  y = quantile(q_i_med,[1/3, 2/3])


figure
hold on
%Repeat this 4 times for each analysis
%Normalized Flow
[idx_new,C,n_clust]=kmediods_best(q_avg);
cluster_subplot_old(idx_new,n_clust,C,q_avg,2,2,1,'Day of Week','Normalized Flow','Normalized Water Consumption',0)
idx_q=idx_new;
%Normalized Energy
[idx_new,C,n_clust]=kmediods_best(e_avg);
cluster_subplot_old(idx_new,n_clust,C,e_avg,2,2,3,'Day of Week','Normalized Energy','Normalized Energy Consumption',0)
idx_e=idx_new;
%Water Intensity
[idx_new,C,n_clust]=kmediods_best(q_avg_qi);
cluster_subplot_old(idx_new,n_clust,C,q_avg_qi,2,2,2,'Day of Week','Water Intensity (gal/sqft)','Water Intensity',0)
idx_qi=idx_new;
%Energy Intensity
[idx_new,C,n_clust]=kmediods_best(q_avg_ei);
cluster_subplot_old(idx_new,n_clust,C,q_avg_ei,2,2,4,'Day of Week','Energy Intensity (kWh/sqft)','Energy Intensity',0)
idx_ei=idx_new;


print(['4 Metric Cluster.jpg'],'-dpdf','-bestfit')

%%

%Put cluster results into table
results=table(builds,idx_q,idx_e,idx_qi,idx_ei);

%% Dissimilarity matrix

%Create matrix of all the idx vectors so it's easier to work with
idx=table(idx_q, idx_e, idx_qi, idx_ei);
idx.Properties.RowNames=builds;

s = RandStream('mlfg6331_64');

%loop for each building 
for k=1:73
disp(k)
%do random sample 100 times
for z=1:5000 
    
%Select reference clusters
for i=1:3
clust_ref(i)=randsample(s,find(idx_q==i),1);
end

%Find which cluster it is in for each of the other clusters to make a
%matrix
for i=1:4
dis_mat(i,:)=idx{clust_ref,i};
end

%Score for one building

%Take the first building as a test
build_inds=idx{k,:}';

%Determine which buildings are in the same matrix as the references
for i=1:3
dis_mat_build(:,i)=build_inds==dis_mat(:,i);
end

mat_sum=sum(dis_mat_build);

dis_mat_tot(z,:)=mat_sum==max(mat_sum);

end
dis_mat_all(k,1:3)=sum(dis_mat_tot);
end
%%
figure
plot(dis_mat_all')

%figure
[idx_new,C,n_clust]=kmediods_best(dis_mat_all);
n_clust=9;
[idx_new,C]=kmedoids(dis_mat_all,n_clust);
%[idx_new,C,n_clust]=kmeans_best(dis_mat_all);
%cluster_subplot(idx_new,n_clust,C,dis_mat_all,1,1,1,'Original Cluster','N','Dissimilarity Matrix',0)

figure
subplot(1,2,1)
plot(dis_mat_all')
cluster_subplot_old(idx_new,n_clust,C,dis_mat_all,1,2,2,'Original Cluster','Number of Simulations','Dissimilarity Matrix',0);






table = Table.open("C:\\Users\\JABER003\\OneDrive - WageningenUR\\Paper\\Heterogeneity\\Data\\EggYolk\\Droplets_coardinate.csv");

n = Table.size
xmin = newArray(n);
xmax = newArray(n);
ymin = newArray(n);
ymax = newArray(n);
x_f = newArray(n);
y_f = newArray(n);
for(i = 0; i < n; i++)
{	correction = 26.04;
	xmin[i] = getResult("xmin", i)*correction; //you need to modify this numbers
	xmax[i] = getResult("xmax", i)*correction;
	ymin[i] = getResult("ymin", i)*correction;
	ymax[i] = getResult("ymax", i)*correction;
	
}
print(n);

run("Import results", "detectmeasurementprotocol=false filepath=[C:\\Users\\JABER003\\OneDrive - WageningenUR\\Paper\\Heterogeneity\\Data\\EggYolk\\egg_yolk_red5_1\\Data.csv] fileformat=[CSV (comma separated)] livepreview=true rawimagestack= startingframe=1 append=false");

for(i = 0; i < n; i++)
{
	run("Show results table", "action=filter formula=[(x > "+ xmin[i] +" & x <"+ xmax[i] +" & y > "+ ymin[i] +" & y < "+ ymax[i] +")]");
	run("Export results", "floatprecision=5 filepath=[C:\\Users\\JABER003\\OneDrive - WageningenUR\\Paper\\Heterogeneity\\Data\\EggYolk\\Red_channel\\Droplet_"+i + 1+".csv] fileformat=[CSV (comma separated)] intensity=false offset=false saveprotocol=false x=true sigma2=false y=true sigma1=true z=false bkgstd=false id=false frame=false");
	run("Show results table", "action=reset");	
}


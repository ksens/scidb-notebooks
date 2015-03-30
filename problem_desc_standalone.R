# Build an array with three attributes per cell
m = scidb("build(<vid:int32>[i=0:100000,1000,0], random()%10000)")
m = bind(m, "x", "bool(random()%2)")
m = bind(m, "y", "bool(random()%2)")
m = scidbeval(m)

# This aggregation fails due to memory shortage
A <- aggregate(m, by=list("vid", "x", "y"), FUN="count(*)", eval=TRUE)

# This aggregation succeeds
B <- aggregate(m, by=list("vid"), FUN="count(*)", eval=TRUE)



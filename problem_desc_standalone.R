# Build an array with three attributes per cell
m = scidb("build(<vid:int32>[i=0:100000,1000,0], random()%10000)")
m = bind(m, "x", "bool(random()%2)")
m = bind(m, "y", "bool(random()%2)")
m = scidbeval(m)

# This aggregation fails due to memory shortage
A <- aggregate(m, by=list("vid", "x", "y"), FUN="count(*)", eval=TRUE)

# This aggregation succeeds
B <- aggregate(m, by=list("vid"), FUN="count(*)", eval=TRUE)

## ---- Rewrite with int64 types ----------
# Build an array with three attributes per cell
m = scidb("build(<vid:int64>[i=0:100000,1000,0], random()%10000)")
m = bind(m, "x", "int64(random()%2)")
m = bind(m, "y", "int64(random()%2)")
m = scidbeval(m)

m1 <- bind(m,"one",1)

# fails with a syntax error; works after updating the scidb-R client
# from github
A <- redimension(m1,
                 dim=c("vid","x","y"),
                 FUN="sum(one) as sum")




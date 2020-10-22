test_that("running without error with defaults",{
expect_error(comm.simul(x=seq(0.1,0.9,0.08),S=100,J=200),regexp=NA)
})

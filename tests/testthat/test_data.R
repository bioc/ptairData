testthat::context("Test data.")

test_data <- function() {
  directory <-  system.file("extdata",  package = "ptairData")
  testthat::expect_equal( list.files(directory) , c("exhaledAir"  ,
                                                    "mycobacteria","simulation"))
  testthat::expect_equal( list.files(file.path(directory,"exhaledAir")) ,
                          c("ind1"  , "ind2"))
  testthat::expect_equal( list.files(file.path(directory,"exhaledAir/ind1")) ,
                          c("ind1-1.h5"  , "ind1-2.h5","ind1-3.h5"))
  testthat::expect_equal( list.files(file.path(directory,"mycobacteria")),
                          c("Control" , "Specie-a", "Specie-b"))
  testthat::expect_equal( list.files(file.path(directory,"mycobacteria"),
                                     recursive = TRUE), c("Control/Control1.h5" ,
                                                          "Control/Control2.h5" ,
                                                          "Specie-a/Specie-a1.h5" ,
                                                          "Specie-a/Specie-a2.h5" ,
                                                          "Specie-b/specie-b1.h5",
                                                          "Specie-b/specie-b2.h5"))

}

testthat::test_that("All data are presents.", test_data())

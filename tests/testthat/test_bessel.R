test_that("BesselK_real_line", {
  z <- seq(-3, 3, by = 0.35)
  bessel_test <- bessel_k(z, 2)
  Bessel_test <- BesselK(z, 2)
  expect_equal(bessel_test, Bessel_test)
})

test_that("BesselK_positive", {
  z <- seq(0.35, 10, by = 0.01)
  bessel_test <- bessel_k(z, 2)
  Bessel_test <- BesselK(z, 2)
  expect_equal(bessel_test, Bessel_test)
})

test_that("BesselK_complex", {
  z <- seq(-10, 10, by = 0.35) + 1i
  bessel_test <- bessel_k(z, 2)
  Bessel_test <- BesselK(z, 2)
  expect_equal(bessel_test, Bessel_test)
})

test_that("BesselK_complex_scaled", {
  z <- seq(-10, 10, by = 0.35) + 1i
  bessel_test <- bessel_k(z, 2, expon_scaled = TRUE)
  Bessel_test <- BesselK(z, 2, expon.scaled = TRUE)
  expect_equal(bessel_test, Bessel_test)
})


test_that("BesselI_real_line", {
  z <- seq(-3, 3, by = 0.35)
  bessel_test <- bessel_i(z, 2)
  Bessel_test <- BesselI(z, 2)
  expect_equal(bessel_test, Bessel_test)
})

test_that("BesselI_positive", {
  z <- seq(0.35, 10, by = 0.01)
  bessel_test <- bessel_i(z, 2)
  Bessel_test <- BesselI(z, 2)
  expect_equal(bessel_test, Bessel_test)
})

test_that("BesselI_complex", {
  z <- seq(-10, 10, by = 0.35) + 1i
  bessel_test <- bessel_i(z, 2)
  Bessel_test <- BesselI(z, 2)
  expect_equal(bessel_test, Bessel_test)
})

test_that("BesselI_complex_scaled", {
  z <- seq(-10, 10, by = 0.35) + 1i
  bessel_test <- bessel_i(z, 2, expon_scaled = TRUE)
  Bessel_test <- BesselI(z, 2, expon.scaled = TRUE)
  expect_equal(bessel_test, Bessel_test)
})



test_that("BesselJ_real_line", {
  z <- seq(-3, 3, by = 0.35)
  bessel_test <- bessel_j(z, 2)
  Bessel_test <- BesselJ(z, 2)
  expect_equal(bessel_test, Bessel_test)
})

test_that("BesselJ_positive", {
  z <- seq(0.35, 10, by = 0.01)
  bessel_test <- bessel_j(z, 2)
  Bessel_test <- BesselJ(z, 2)
  expect_equal(bessel_test, Bessel_test)
})

test_that("BesselJ_complex", {
  z <- seq(-10, 10, by = 0.35) + 1i
  bessel_test <- bessel_j(z, 2)
  Bessel_test <- BesselJ(z, 2)
  expect_equal(bessel_test, Bessel_test)
})

test_that("BesselJ_complex_scaled", {
  z <- seq(-10, 10, by = 0.35) + 1i
  bessel_test <- bessel_j(z, 2, expon_scaled = TRUE)
  Bessel_test <- BesselJ(z, 2, expon.scaled = TRUE)
  expect_equal(bessel_test, Bessel_test)
})


test_that("BesselY_real_line", {
  z <- seq(-3, 3, by = 0.35)
  bessel_test <- bessel_y(z, 2)
  Bessel_test <- BesselY(z, 2)
  expect_equal(bessel_test, Bessel_test)
})

test_that("BesselY_positive", {
  z <- seq(0.35, 10, by = 0.01)
  bessel_test <- bessel_y(z, 2)
  Bessel_test <- BesselY(z, 2)
  expect_equal(bessel_test, Bessel_test)
})

test_that("BesselY_complex", {
  z <- seq(-10, 10, by = 0.35) + 1i
  bessel_test <- bessel_y(z, 2)
  Bessel_test <- BesselY(z, 2)
  expect_equal(bessel_test, Bessel_test)
})

test_that("BesselY_complex_scaled", {
  z <- seq(-10, 10, by = 0.35) + 1i
  bessel_test <- bessel_y(z, 2, expon_scaled = TRUE)
  Bessel_test <- BesselY(z, 2, expon.scaled = TRUE)
  expect_equal(bessel_test, Bessel_test)
})


test_that("BesselH_real_line", {
  z <- seq(-3, 3, by = 0.35)
  bessel_test <- bessel_h(1, z, 2)
  Bessel_test <- BesselH(1, z, 2)
  expect_equal(bessel_test, Bessel_test)
})

test_that("BesselH_positive", {
  z <- seq(0.35, 10, by = 0.01)
  bessel_test <- bessel_h(1, z, 2)
  Bessel_test <- BesselH(1, z, 2)
  expect_equal(bessel_test, Bessel_test)
})

test_that("BesselH_complex", {
  z <- seq(-10, 10, by = 0.35) + 1i
  bessel_test <- bessel_h(1, z, 2)
  Bessel_test <- BesselH(1, z, 2)
  expect_equal(bessel_test, Bessel_test)
})

test_that("BesselH_complex_scaled", {
  z <- seq(-10, 10, by = 0.35) + 1i
  bessel_test <- bessel_h(1, z, 2, expon_scaled = TRUE)
  Bessel_test <- BesselH(1, z, 2, expon.scaled = TRUE)
  expect_equal(bessel_test, Bessel_test)
})


test_that("BesselH_real_line_2", {
  z <- seq(-3, 3, by = 0.35)
  bessel_test <- bessel_h(2, z, 2)
  Bessel_test <- BesselH(2, z, 2)
  expect_equal(bessel_test, Bessel_test)
})

test_that("BesselH_positive_2", {
  z <- seq(0.35, 10, by = 0.01)
  bessel_test <- bessel_h(2, z, 2)
  Bessel_test <- BesselH(2, z, 2)
  expect_equal(bessel_test, Bessel_test)
})

test_that("BesselH_complex_2", {
  z <- seq(-10, 10, by = 0.35) + 1i
  bessel_test <- bessel_h(2, z, 2)
  Bessel_test <- BesselH(2, z, 2)
  expect_equal(bessel_test, Bessel_test)
})

test_that("BesselH_complex_scaled_2", {
  z <- seq(-10, 10, by = 0.35) + 1i
  bessel_test <- bessel_h(2, z, 2, expon_scaled = TRUE)
  Bessel_test <- BesselH(2, z, 2, expon.scaled = TRUE)
  expect_equal(bessel_test, Bessel_test)
})




test_that("AiryA_real_0", {
  z <- seq(-3, 3, by = 0.35)
  airy_test <- airy_a(z, 0)
  Airy_test <- AiryA(z, 0)
  expect_equal(airy_test, Airy_test)
})


test_that("AiryA_complex_0", {
  z <- seq(-3, 3, by = 0.35) + 1i
  airy_test <- airy_a(z, 0)
  Airy_test <- AiryA(z, 0)
  expect_equal(airy_test, Airy_test)
})


test_that("AiryA_real_1", {
  z <- seq(-3, 3, by = 0.35)
  airy_test <- airy_a(z, 1)
  Airy_test <- AiryA(z, 1)
  expect_equal(airy_test, Airy_test)
})


test_that("AiryA_complex_1", {
  z <- seq(-3, 3, by = 0.35) + 1i
  airy_test <- airy_a(z, 1)
  Airy_test <- AiryA(z, 1)
  expect_equal(airy_test, Airy_test)
})


test_that("AiryB_real_0", {
  z <- seq(-3, 3, by = 0.35)
  airy_test <- airy_b(z, 0)
  Airy_test <- AiryB(z, 0)
  expect_equal(airy_test, Airy_test)
})


test_that("AiryB_complex_0", {
  z <- seq(-3, 3, by = 0.35) + 1i
  airy_test <- airy_b(z, 0)
  Airy_test <- AiryB(z, 0)
  expect_equal(airy_test, Airy_test)
})


test_that("AiryB_real_1", {
  z <- seq(-3, 3, by = 0.35)
  airy_test <- airy_b(z, 1)
  Airy_test <- AiryB(z, 1)
  expect_equal(airy_test, Airy_test)
})


test_that("AiryB_complex_1", {
  z <- seq(-3, 3, by = 0.35) + 1i
  airy_test <- airy_b(z, 1)
  Airy_test <- AiryB(z, 1)
  expect_equal(airy_test, Airy_test)
})




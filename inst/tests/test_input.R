context("tests on inputs")

test_that("processPheno works", {
      expect_error(processPheno("SRP00000"), "Invalid 'project' argument. There's no such 'project' in the recount_url data.frame.")

      project = 'SRP063581'
      pheno = processPheno(project)
      expect_equal(dim(pheno)[1], 3)
      expect_equal(dim(pheno)[2], 24)
      expect_equal(pheno$project[1], 'SRP063581')
      expect_equal(pheno$avg_read_length, pheno$rls)

      pheno = processPheno("ERP001942")
      expect_equal(dim(pheno)[1], 663)
      expect_equal(dim(pheno)[2], 24)
      expect_equal(pheno$project[1], 'ERP001942')
      expect_equal(pheno$avg_read_length[1], 150)
      expect_equal(pheno$rls[1], 75)

      input = data.frame(project="Test", run=1:5, bigwig_path=paste0("path", 1:5), rls = 50)
      expect_error(processPheno(input))

      input = data.frame(project="Test", run=1:5, bigwig_path=paste0("path", 1:5), rls = "50", paired_end=T)
      expect_error(processPheno(input), "rls should be numeric")

      input = data.frame(project="Test", run=1:5, bigwig_path=paste0("path", 1:5), rls = 50, paired_end=10)
      expect_error(processPheno(input), 'paired_end must be T/F')

      input = data.frame(project="Test", run=1:5, bigwig_path=paste0("path", 1:5), rls = 50, paired_end=T)
      pheno = processPheno(input)
      expect_equal(dim(pheno)[1], 5)
      expect_equal(dim(pheno)[2], 6)
})

test_that("getJxCounts works", {
      pheno = processPheno("SRP002915")
      expect_equal(getJxCounts("SRP002915", pheno), NULL)

      pheno = processPheno("DRP000366")
      counts_jx = getJxCounts("DRP000366", pheno)
      expect_equal(counts_jx[1:6], c(5, 2, 6, 2, 1, 16))
})

test_that("getExCounts works", {
      pheno = processPheno("DRP000366")
      counts_ex = getExCounts(pheno)
      expect_equal(dim(counts_ex)[1], 2529474)
      expect_equal(dim(counts_ex)[2], 1)
      expect_equal(sum(counts_ex[,1]), 18872281)
})

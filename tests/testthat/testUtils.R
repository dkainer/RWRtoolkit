context("Utils Tests")
library(RWRtoolkit)
library(vctrs)
library(mockery)
load('../testMultiplex/unitTestMultiplex.Rdata')

describe('load_geneset', {
    it('loads genes with sharps within the names', {
        # TODO: 
        # Add fill*, trunctes the remainder of the line
        # commend character paramter - metioned by @Izaak comment.char='#'
    })


    it('returns null if geneset path is null', {       
        output <- RWRtoolkit::load_geneset(NULL)

        expect_equal(output, NULL)
    })

    it('fails to load the geneset due to bad file path', {
        filePath <- '../testGenesets/bad/file/path.tsv'

        expectedErrorMessage <- paste("ERROR: geneset file does not exist:", filePath)
        expect_error(RWRtoolkit::load_geneset(filePath), expectedErrorMessage)
    })

    it('loads bad geneset', {
        filePath <- '../testGenesets/testGeneset_Bad_oneColumn.tsv'
        expectedErrorMessage <- "Your geneset file is incorrectly formatted. Please see documentation."

        expect_error(RWRtoolkit::load_geneset(filePath), expectedErrorMessage)
    })

    it('loads incorrectly formatted geneset', {
        filePath <- '../testGenesets/testGeneset1.csv'
        expectedErrorMessage <- "Your geneset file is incorrectly formatted. Please see documentation."

        expect_error(RWRtoolkit::load_geneset(filePath), expectedErrorMessage)
    })

    it('loads unfiltered geneset with no multiplex network', {
        filePath <- '../testGenesets/testGeneset1.tsv'
        setid <- c('setA', 'setA')
        gene <- c('1','2')
        weight <- c(1, 1)
        expectedGeneset <- data.frame(setid, gene, weight)
        expectedExtras <- NULL        
        expectedOutput <- list("geneset"=expectedGeneset, "extras"=expectedExtras)

        output <- RWRtoolkit::load_geneset(filePath)
        
        expect_equal(output, expectedOutput)
    })

    it('loads weighted unfiltered geneset  with no multiplex network', {
        filePath <- '../testGenesets/testGeneset2.tsv'
        setid <- c('setB', 'setB')
        gene <- c('2','3')
        weight <- c(0.4, 0.6)
        expectedGeneset <- data.frame(setid, gene, weight)
        expectedExtras <- NULL
        expectedOutput <- list("geneset"=expectedGeneset, "extras"=expectedExtras)

        output <- RWRtoolkit::load_geneset(filePath)

        expect_equal(output, expectedOutput)
    })

    it('loads un-weighted unfiltered geneset with 3 columns with no multiplex network', {
        ## Has a third column, but words, not weights. 
        filePath <- '../testGenesets/testGeneset3.tsv'
        setid <- c('setC', 'setC')
        gene <- c('4','5')
        weight <- c(1, 1)
        expectedGeneset <- data.frame(setid, gene, weight)
        expectedExtras <- NULL
        expectedOutput <- list("geneset"=expectedGeneset, "extras"=expectedExtras)

        output <- RWRtoolkit::load_geneset(filePath)

        expect_equal(output, expectedOutput)
    })
    
    it('loads filtered geneset filtered by the multiplex pool of nodes', {  
        # Even though testGeneset2 has two genes: (2, 3) - the pool of nodes has nodes 0, 1, 2. 
        # The intersection of the two is specifically gene 2, which is the only output we expect to see. 
        Pool_of_Nodes <- c('0','1','2') 
        nw.mpo <- data.frame(Pool_of_Nodes)
        filePath <- '../testGenesets/testGeneset2.tsv'
        setid <- c('setB')
        gene <- c('2')
        weight <- c(0.4)
        expectedGeneset <- data.frame(setid, gene, weight)
        gene <- c('3')
        weight <- c(0.6)
        expectedExtras <- data.frame(setid, gene, weight)
        expectedOutput <- list("geneset"=expectedGeneset, "extras"=expectedExtras)

        output <- expect_warning(RWRtoolkit::load_geneset(filePath, nw.mpo))
        expect_equal(output, expectedOutput)
    })

    it('loads filtered geneset filtered by the multiplex pool of nodes, all filtered', {
        # Even though testGeneset2 has two genes: (2, 3) - the pool of nodes has nodes 0, 1, 2. 
        # The intersection of the two is specifically gene 2, which is the only output we expect to see. 
        Pool_of_Nodes <- c('0','1') 
        nw.mpo <- data.frame(Pool_of_Nodes)
        filePath <- '../testGenesets/testGeneset2.tsv'
        expectedGeneset <- data.frame(setid=character(), gene=character(), weight=double())
        setid <- c('setB','setB')
        gene <- c('2','3')
        weight <- c(0.4,0.6)
        expectedExtras <- data.frame(setid, gene,weight)
        expectedOutput <- list("geneset"=expectedGeneset, "extras"=expectedExtras)

        output <- expect_warning(RWRtoolkit::load_geneset(filePath, nw.mpo))

        expect_equal(output, expectedOutput)
    })
})

describe('get_or_set_tau', {
    it('Returns numeric value if tau is numeric', {
        # Where tau values for layer one = 1.25 and layer 2 = 0.75

        optTau <- c(1.25, 0.75)
        expectedOutput <- optTau

        output <- RWRtoolkit::get_or_set_tau(nw.mpo, optTau)

        expect_equal(output, expectedOutput)
    })

    it('Returns numeric values if tau is string list obtained via command line', {
        # Where tau values for layer one = 1.25 and layer 2 = 0.75
        optTau <- '1.25,0.75'
        expectedOutput <- c(1.25, 0.75)

        output <- RWRtoolkit::get_or_set_tau(nw.mpo, optTau)

        expect_equal(output, expectedOutput)
    })

    it('Returns one for each layer with warning because of incorrect sum of tau vals', {
        # Where tau values for layer one = 0.25 and layer 2 = 0.75, don't sum to 2
        optTau <- c(0.25, 0.75)
        expectedOutput <- c(1,1)

        output <- expect_warning(RWRtoolkit::get_or_set_tau(nw.mpo, optTau))

        expect_equal(output, expectedOutput)
    })

    it('Returns one for each layer with warning because of incorrect list size', {
        # Where tau values for layer one = 0.25 and layer 2 = 0.75, don't sum to 2
        optTau <- '0.75'
        expectedOutput <- c(1,1)

        output <- expect_warning(RWRtoolkit::get_or_set_tau(nw.mpo, optTau))

        expect_equal(output, expectedOutput)
    })    
})

describe('chunk', {
    it('takes in a dataframe, and separates it by N values (number of folds for kfold cross val)', {        
        geneList <- c('a','b','c','d','e','f')
        nFolds  <-  3
        expectedOutputLength  <-  3

        output <- RWRtoolkit::chunk(geneList, nFolds)

        expect_equal(length(names(output)), expectedOutputLength)
        expect_equal(output['1'][[1]], c('a','b'))
        expect_equal(output['2'][[1]], c('c','d'))
        expect_equal(output['3'][[1]], c('e','f'))
    })

    it('Warns users if gene set length is smaller than num folds', {
        geneList <- c('a','b','c')
        nFolds  <-  4
        expectedOutputLength  <-  3

        output <- expect_warning(RWRtoolkit::chunk(geneList, nFolds))

        expect_equal(length(names(output)), expectedOutputLength)
    })
})

describe('write_table', {
    it('Fails to write table if no table exists', {
        table <-  NA
        path <-  './tableOutput'

        expect_warning(RWRtoolkit::write_table(table, path))
    })

    it('Writes table to existing path', {
        fauxListOfNodes <-  c('A', 'B', 'C', 'D')        
        path <-  './someFakeFile.txt'

        stubWriteTable <-  mock()
        stub(RWRtoolkit::write_table, 'write.table', stubWriteTable)
        RWRtoolkit::write_table(fauxListOfNodes, path)

        expect_called(stubWriteTable, 1)
        expect_args(stubWriteTable, 1, fauxListOfNodes, path, '\t', F, T, F)
    })

    it('Writes table to nonexisting path', {
        fauxListOfNodes <-  c('A', 'B', 'C', 'D')        
        nonExistentPath <-  './someOtherPath/recursiveDepth/output.txt'
        expectedDirectory <- 'recursiveDepth'

        stubWriteTable <-  mock()
        stub(RWRtoolkit::write_table, 'write.table', stubWriteTable)
        RWRtoolkit::write_table(fauxListOfNodes, nonExistentPath)
        checkDir <-  system('ls ./someOtherPath/', intern=T)

        expect_equal(checkDir[1], expectedDirectory)
        expect_called(stubWriteTable, 1)
        expect_args(stubWriteTable, 1, fauxListOfNodes, nonExistentPath, '\t', F, T, F)
    })
})

describe('get_file_path', {
    it('appends default extension to a filename with no specified out directory', {
        baseFileName <- 'someFile'
        expectedFileName <- 'someFile.tsv'

        outputFileName <- RWRtoolkit::get_file_path(baseFileName)

        expect_equal(outputFileName, expectedFileName)
    })


    it('appends extension to a filename with no specified out directory', {
        baseFileName <- 'someFile'
        expectedFileName <- 'someFile.csv'
        fileExtension <- '.csv'

        outputFileName <- RWRtoolkit::get_file_path(baseFileName, ext=fileExtension)

        expect_equal(outputFileName, expectedFileName)
    })
    it('appends extension to a filename with specified out directory', {
        baseFileName <- 'someFile'
        outDirectoryPath <- '/some/path/to/file'
        expectedFileName <- '/some/path/to/file/someFile.tsv'

        outputFileName <- RWRtoolkit::get_file_path(baseFileName, outdir=outDirectoryPath)

        expect_equal(outputFileName, expectedFileName)
    })
})

describe('dump_layers', {
    it('extracts specific layers and writes them to a table', {
        load('../testMultiplex/network.Rdata')
        fakeFilePath <- 'someFilePath/network.tsv'
        expectedLayer1Name <- 'm1'
        expectedLayer2Name <- 'm2'
        expectedNumberOfCalls <- 2
        from <- c('0', '0', '1')
        to <- c('1', '2', '2')
        weight <- c(0.5, 1.0, 0.5)
        type <- c('m1','m1','m1')
        expectedLayer1DF <- data.frame(from, to, weight, type)
        from <- c('0', '0', '1', '2')
        to <- c('2', '3', '2', '3')
        weight <- c(1, 1, 1, 1)
        type <- c('m2', 'm2', 'm2', 'm2')
        expectedLayer2DF <- data.frame(from, to, weight, type)

        stubGetFilePath <- mock(fakeFilePath, fakeFilePath)
        stubWriteTable <-  mock()
        stub(RWRtoolkit::dump_layers, 'get_file_path', stubGetFilePath)
        stub(RWRtoolkit::dump_layers, 'write_table', stubWriteTable)
        RWRtoolkit::dump_layers(nw.mpo)

        expect_called(stubGetFilePath, expectedNumberOfCalls)
        expect_args(stubGetFilePath, 1, expectedLayer1Name, NULL)
        expect_args(stubGetFilePath, 2, expectedLayer2Name, NULL)
        expect_called(stubWriteTable, 2)
        expect_args(stubWriteTable, 1, expectedLayer1DF, fakeFilePath)
        expect_args(stubWriteTable, 2, expectedLayer2DF, fakeFilePath)
    })
    
    it('writes nothing to a table if no layers exist in object', {
        load('../testMultiplex/network.Rdata')
        expectedNumberOfCalls <- 0
        fakeFilePath <- 'someFilePath/network.tsv'
        nw.mpo$m1 <- NULL # Remove layers from network. 
        nw.mpo$m2 <- NULL

        stubGetFilePath <- mock(fakeFilePath, fakeFilePath)
        stubWriteTable <-  mock()
        stub(RWRtoolkit::dump_layers, 'get_file_path', stubGetFilePath)
        stub(RWRtoolkit::dump_layers, 'write_table', stubWriteTable)
        RWRtoolkit::dump_layers(nw.mpo)

        expect_called(stubGetFilePath, expectedNumberOfCalls)
        expect_called(stubWriteTable, expectedNumberOfCalls)
    })
})

describe('dump_nodes', {
    it('throws an error if a multiplexMH object is not passed in as mpo', {
        nw.mpo <- NULL
        expectedNumberOfCalls <- 0

        stubWriteTable <-  mock()
        stub(RWRtoolkit::dump_nodes, 'write_table', stubWriteTable)
        expect_error(RWRtoolkit::dump_nodes(nw.mpo))

        expect_called(stubWriteTable, expectedNumberOfCalls)
    })
    it('extracts nodes from the multiplexMH object and saves to file', {
        load('../testMultiplex/network.Rdata')
        mpo <- nw.mpo1 # loaded from the network.Rdata
        mpoNodes <- c('0','1','2')
        expectedDataFrame <- data.frame(mpoNodes)
        expectedNumberOfCalls <- 1
        fakeFilePath<- './someFile/path/file.txt'

        stubWriteTable <-  mock()
        stubGetFilePath <- mock(fakeFilePath)
        stub(RWRtoolkit::dump_nodes, 'write_table', stubWriteTable)
        stub(RWRtoolkit::dump_nodes, 'get_file_path', stubGetFilePath)
        RWRtoolkit::dump_nodes(nw.mpo1)
        
        expect_called(stubGetFilePath, expectedNumberOfCalls)
        expect_called(stubWriteTable, expectedNumberOfCalls)
        expect_args(stubWriteTable, 1,  expectedDataFrame, fakeFilePath)
    })
})


# removes test files so as not to clutter everything up
teardown({
    system('rm -rf ./someOtherPath') # removes path to written table via write_table
    system('rm ./tableOutput')
}, env=parent.frame())

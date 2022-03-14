context("Random Walk Restart Tests")
library(RWRtoolkit)
library(RandomWalkRestartMH)
library(vctrs)
library(igraph)
library(mockery)

# Homogenous Network tibble
nwTibble <-tibble::tibble(
    'nwfile'=c('../testNetworks/m1.txt', '../testNetworks/m2.txt'), 
    'nwname'=c('m1','m2'), 
    'nwgroup'=c(1, 1)) 

nwTibbleBadPath <-tibble::tibble(
    'nwfile'=c('../testNetworks/iDontExist_m1.txt', '../whereever/m2.txt'), 
    'nwname'=c('m1','m2'), 
    'nwgroup'=c(1, 1)) 

nw.groups <- list_of(nwTibble)

# Heterogeneous Network Tibble
nwTibble_Het1 <-tibble::tibble(
    'nwfile'=c('../testNetworks/m1.txt'), 
    'nwname'=c('m1'), 
    'nwgroup'=c( 1 )) 
nwTibble_Het2 <-tibble::tibble(
    'nwfile'=c( '../testNetworks/n1.txt'), 
    'nwname'=c('n1'), 
    'nwgroup'=c(2)) 
nwTibble_Het3 <-tibble::tibble(
    'nwfile'=c( '../testNetworks/i1.txt'), 
    'nwname'=c('i1' ), 
    'nwgroup'=c( 3 )) 


# generates expected matrix for test data from graphs m1 and m2
generateExpectedSupraAdjacency <- function(delta){
    oneMinusDelta <- 1 - delta
    # Because of how RandomWalkRestartMH Calculates the edge weights
    # our graph, m1, originally has weights 1 and 2. Those get normalized
    # my the create multiplex to then have the weights 0.5 and 1, respectively.
    w1 <- 1 * oneMinusDelta
    w1Half <- 0.5 * oneMinusDelta

    colNames <- c("0_1", "1_1",  "2_1",  "3_1",  "0_2",  "1_2",  "2_2",  "3_2")
    listMat <- c(
        c(0     ,w1Half ,w1     ,0     ,delta ,0     ,0     ,0    ),
        c(w1Half,0      ,w1Half ,0     ,0     ,delta ,0     ,0    ),
        c(w1    ,w1Half ,0      ,0     ,0     ,0     ,delta ,0    ),
        c(0     ,0      ,0      ,0     ,0     ,0     ,0     ,delta),
        c(delta ,0      ,0      ,0     ,0     ,0     ,w1    ,w1   ),
        c(0     ,delta  ,0      ,0     ,0     ,0     ,w1    ,0    ),
        c(0     ,0      ,delta  ,0     ,w1    ,w1    ,0     ,w1   ),
        c(0     ,0      ,0      ,delta ,w1    ,0     ,w1    ,0    )
    )
    mat <- matrix(listMat, nrow=8, ncol=8, dimnames=list(colNames, colNames))
    #Normalize the matrix by column
    normalizedMatrix <- mat %*% diag(1/colSums(mat))
    #non normalized matrix retains the original delta as the intra-layer weight
    expected_nonNormalizedMatrix <- as(mat, 'dgCMatrix')
    #non matrix are normalized along with the [0,1] normalized edge weights in each layer. 
    expected_normalizedMatrix <- as(normalizedMatrix, 'dgCMatrix')
    colnames(expected_normalizedMatrix) <- colNames
    
    return(c(expected_nonNormalizedMatrix, expected_normalizedMatrix))
}

runTestForDifferingGraphData <- function(nw.groups, delta, outputFileName, verbose) {
    supraAdjMatrices <- generateExpectedSupraAdjacency(delta)
    expected_nonNormalizedMatrix <- supraAdjMatrices[[1]]
    expected_normalizedMatrix <- supraAdjMatrices[[2]]

    invisible(RWRtoolkit::make_homogenous_network(nw.groups, delta, outputFileName, verbose))
    load(outputFileName)

    expect_equal(nw.adjnorm, expected_normalizedMatrix)
    expect_equal(nw.adj, expected_nonNormalizedMatrix)
}

describe('make_multiplex', {
    it('throws an error when fed an flist tibble with non-existant files', {
        nw.groups <- list_of(nwTibbleBadPath)
        expect_error(RWRtoolkit::make_multiplex(nw.groups[[1]]))
    })

    it('also makes a multiplex', {
        #expected output as multiplex object: 
        # Because of how RandomWalkRestartMH Calculates the edge weights
        # our graph, m1, originally has weights 1 and 2. Those get normalized
        # my the create multiplex to then have the weights 0.5 and 1, respectively.
        graph1 <- make_graph(c("0","1","0","2","1","2"), directed=FALSE)
        E(graph1)$weight <- c(0.5, 1.0, 0.5)
        E(graph1)$type <- c("m1", "m1", "m1")
        # Add vertex 3 because RandomWalkRestartMH 
        # ensures all layers share common vertices
        graph1 <- add_vertices(graph1, 1, name=c("3")) 
        graph2 <- make_graph(c("0","2","0","3","1","2", "2", "3"), directed=FALSE)
        E(graph2)$weight <- c(1, 1, 1, 1)
        E(graph2)$type <- c("m2", "m2", "m2", "m2")
        expectedNumLayers <- 2
        expectedNumNodes <- 4

        invisible(capture.output(output <- RWRtoolkit::make_multiplex(nw.groups[[1]])))

        expect_equal(output$Number_of_Layers, expectedNumLayers)
        expect_equal(output$Number_of_Nodes_Multiplex, expectedNumNodes)
        # Issues with sorted nature of vertices and edges.   
        # Check node equality for each layer: 
        expect_true(all(sort(V(output$m1)) == sort(V(graph1))))
        expect_true(all(sort(V(output$m2)) == sort(V(graph2))))
        # Check edge equality for each layer: 
        expect_true(all(sort(E(output$m1)) == sort(E(graph1))))
        expect_true(all(sort(E(output$m2)) == sort(E(graph2))))
        # Check edge weight equality for each layer: 
        expect_equal(E(output$m1)$weight, E(graph1)$weight)
        expect_equal(E(output$m2)$weight, E(graph2)$weight)
        # Check Layer descriptions for each layer: 
        expect_equal(E(output$m1)$type, E(graph1)$type)
        expect_equal(E(output$m2)$type, E(graph2)$type)
   })
})

describe('make_homogenous_network', {
    it('takes in an nw.groups object and saves multiplex network data to file', {
        # delta is the non-normalized edge weight of the connecting edges 
        delta <- 0.5
        outputFileName <- 'testthatOutput.txt'
        verbose <- FALSE

        invisible(RWRtoolkit::make_homogenous_network(nw.groups, delta, outputFileName, verbose))

        expect_true(outputFileName %in% list.files())
    })
    
    
    it('ensures multiplex save file matches expected data', { 
        delta <- 0.5
        outputFileName <- 'testthatOutputP5.txt'
        verbose <- FALSE
        
        runTestForDifferingGraphData(nw.groups, delta, outputFileName, verbose)
    })
    
    it('ensures multiplex save file matches expected data with 0.7 as delta', {
        delta <- 0.7
        outputFileName <- 'testthatOutputP7.txt'
        verbose <- FALSE
        
        runTestForDifferingGraphData(nw.groups, delta, outputFileName, verbose)
    })

    it('ensures multiplex save file matches expected data with 1.0 as delta', {
        outputFileName <- 'testthatOutput1P0.txt'
        ## By calling 1 as our delta, the network layers m1 and m2 are entirely 
        ## ignored, and only our delta intra-layer edges exist. 
        delta <- 1
        colNames <- c("0_1", "1_1",  "2_1",  "3_1",  "0_2",  "1_2",  "2_2",  "3_2")
        listMat <- c(
            c(0     ,0      ,0      ,0     ,delta ,0     ,0     ,0    ),
            c(0     ,0      ,0      ,0     ,0     ,delta ,0     ,0    ),
            c(0     ,0      ,0      ,0     ,0     ,0     ,delta ,0    ),
            c(0     ,0      ,0      ,0     ,0     ,0     ,0     ,delta),
            c(delta ,0      ,0      ,0     ,0     ,0     ,0     ,0    ),
            c(0     ,delta  ,0      ,0     ,0     ,0     ,0     ,0    ),
            c(0     ,0      ,delta  ,0     ,0     ,0     ,0     ,0    ),
            c(0     ,0      ,0      ,delta ,0     ,0     ,0     ,0    )
        )
        mat <- matrix(listMat, nrow=8, ncol=8, dimnames=list(colNames, colNames))
        #Normalize the matrix by column
        normalizedMatrix <- mat %*% diag(1/colSums(mat))
        #non normalized matrix retains the original delta as the intra-layer weight
        expected_nonNormalizedMatrix <- as(mat, 'dgCMatrix')
        #non matrix are normalized along with the [0,1] normalized edge weights in each layer. 
        expected_normalizedMatrix <- as(normalizedMatrix, 'dgCMatrix')
        colnames(expected_normalizedMatrix) <- colNames

        invisible(RWRtoolkit::make_homogenous_network(nw.groups, delta, outputFileName, verbose))
        load(outputFileName)

        expect_equal(nw.adj, expected_nonNormalizedMatrix)
        expect_equal(nw.adjnorm, expected_normalizedMatrix)
    })

    it('creates multiplexes with tab delimited data', {
        nwTibble <-tibble::tibble(
            'nwfile'=c('../testNetworks/m1_tabDelim.txt', '../testNetworks/m2_tabDelim.txt'), 
            'nwname'=c('m1','m2'), 
            'nwgroup'=c(1, 1)) 

        nw.groups.tabDelimited <- list_of(nwTibble)
        delta <- 0.5
        outputFileName <- 'testthatOutputP5_fromTabDelmited.txt'
        verbose <- FALSE
  
        runTestForDifferingGraphData(nw.groups.tabDelimited, delta, outputFileName, verbose)
    })

    it('creates multiplexes with non-header data', {
        nwTibble <-tibble::tibble(
            'nwfile'=c('../testNetworks/m1_noHeader.txt', '../testNetworks/m2_noHeader.txt'), 
            'nwname'=c('m1','m2'), 
            'nwgroup'=c(1, 1)) 

        nw.groups.noHeader <- list_of(nwTibble)
        delta <- 0.5
        outputFileName <- 'testthatOutputP5_fromNoHeader.txt'
        verbose <- FALSE
  
        runTestForDifferingGraphData(nw.groups.noHeader, delta, outputFileName, verbose)
    })

    it('creates multiplexes with mixed delimited data', {
        nwTibble <-tibble::tibble(
            'nwfile'=c('../testNetworks/m1.txt', '../testNetworks/m2_tabDelim.txt'), 
            'nwname'=c('m1','m2'), 
            'nwgroup'=c(1, 1)) 
        
        nw.groups.mixedDelim <- list_of(nwTibble)
        delta <- 0.5
        outputFileName <- 'testthatOutputP5_fromMixedDelim.txt'
        verbose <- FALSE
  
        runTestForDifferingGraphData(nw.groups.mixedDelim, delta, outputFileName, verbose)
    })
    it('creates multiplexes with mixed header data', {
        nwTibble <-tibble::tibble(
            'nwfile'=c('../testNetworks/m1_noHeader.txt', '../testNetworks/m2.txt'), 
            'nwname'=c('m1','m2'), 
            'nwgroup'=c(1, 1)) 

        nw.groups.mixedHeader <- list_of(nwTibble)
        delta <- 0.5
        outputFileName <- 'testthatOutputP5_fromMixedHeader.txt'
        verbose <- FALSE
  
        runTestForDifferingGraphData(nw.groups.mixedHeader, delta, outputFileName, verbose)
    })
    it('creates multiplexes with mixed header and mixed delimiter data', {
        nwTibble <-tibble::tibble(
            'nwfile'=c('../testNetworks/m1_noHeader.txt', '../testNetworks/m2_tabDelim.txt'), 
            'nwname'=c('m1','m2'), 
            'nwgroup'=c(1, 1)) 

        nw.groups.mixedHeaderDelim <- list_of(nwTibble)
        delta <- 0.5
        outputFileName <- 'testthatOutputP5_fromMixedHeaderDelim.txt'
        verbose <- FALSE
  
        runTestForDifferingGraphData(nw.groups.mixedHeaderDelim, delta, outputFileName, verbose)
    })

    it('throws an error when it fails to create multiplex with a delta value of 0.0', {
        # make_homogenous_network <- function(nw.groups, delta, out, verbose) {
        delta <- 0.0
        oneMinusDelta <- 1 - delta
        outputFileName <- 'testthatOutputP7.txt'
        verbose <- FALSE
        
        expect_error(RWRtoolkit::make_homogenous_network(nw.groups, delta, outputFileName, verbose))
    })


})

describe('make_heterogeneous_multiplex', {
    it('creates a heterogeneous multiplex and saves it to file', {
        ############ QUESTIONS FOR DAVID WITH MAKE MULTIPLEX HET
        
        #         0_1  1_1       2_1 3_1       4_1 5_1 6_1 7_1 8_1
        # 0_1 .         0.25 0.3333333 0.5 .         0.5 .     .   .
        # 1_1 0.1666667 .    0.1666667 .   0.5000000 .   .     .   .
        # 2_1 0.3333333 0.25 .         .   .         .   0.5   .   .
        # 3_1 0.2500000 .    .         .   0.1666667 .   .     .   .
        # 4_1 .         0.50 .         0.5 .         0.5 .     .   1
        # 5_1 0.2500000 .    .         .   0.1666667 .   .     .   .
        # 6_1 .         .    0.5000000 .   .         .   .     1   .
        # 7_1 .         .    .         .   .         .   0.5   .   .
        # 8_1 .         .    .         .   0.1666667 .   .     .   .


        ## Inter layer matrix edge weidhts depend on the labda. 
        ## When given a 0.5, if there existss only one interlayer edge
        ## That edge recieves the weight of lambda. The remainder of the 
        ## weights that are not lambda are then calculated 1-lambda 
        nw.groupsInput <- list_of( nwTibble_Het1,nwTibble_Het2, nwTibble_Het3 )
        delta <- 1

        ## currently, with package formation, delta appears to do nothing? 

        lambda <- 0.6
        out <- "network.Rdata"

        invisible(RWRtoolkit::make_heterogeneous_multiplex(nw.groupsInput, delta, lambda, out))

    })

})

describe('RWR_make_multiplex.R:', {
    it('throws an error if flist is empty', {
        flist <- ''
        
        expect_error(RWRtoolkit::RWR_make_multiplex(flist))
    })

    it('throws an error if flist elements contain bad paths', {
        badFlist <- '../testFlists/testFlist_badPaths.txt'
        
        expect_error(RWRtoolkit::RWR_make_multiplex(badFlist))
    })

    it('takes flist and makes a homogenous multiplex with default parameters', {
        nw.groupsInput <- list_of(nwTibble)
        flistFilePath <- '../testFlists/testFlist.txt'

        makeHomogenousStub = mock()
        makeHeterogenousStub = mock()
        stub(RWRtoolkit::RWR_make_multiplex, 'make_homogenous_network', makeHomogenousStub)
        stub(RWRtoolkit::RWR_make_multiplex, 'make_heterogeneous_multiplex', makeHeterogenousStub)
        
        invisible(RWRtoolkit::RWR_make_multiplex(flistFilePath))

        expect_called(makeHomogenousStub, 1)
        expect_args(makeHomogenousStub, 1, nw.groupsInput, 0.5, 'network.Rdata', FALSE)
        expect_called(makeHeterogenousStub, 0)
    })

    it('takes flist and makes a homogenous multiplex with non-default parameters', {
        nw.groupsInput <- list_of(nwTibble)
        flistFilePath <- '../testFlists/testFlist.txt'

        delta <- 0.9
        verbose <- TRUE
        outputFile <- 'myOutput.txt'
        
        makeHomogenousStub = mock()
        makeHeterogenousStub = mock()
        stub(RWRtoolkit::RWR_make_multiplex, 'make_homogenous_network', makeHomogenousStub)
        stub(RWRtoolkit::RWR_make_multiplex, 'make_heterogeneous_multiplex', makeHeterogenousStub)
        
        invisible(RWRtoolkit::RWR_make_multiplex(flistFilePath, delta=delta, output=outputFile, verbose=verbose))

        expect_called(makeHomogenousStub, 1)
        expect_args(makeHomogenousStub, 1, nw.groupsInput, delta, outputFile, verbose)
        expect_called(makeHeterogenousStub, 0)
    })
    
    it('fails to create a heterogeneous network due to there not being enough files', {
        flistFilePath <- '../testFlists/testFlist_heterogeneous_badGrouping.txt'

        makeHomogenousStub = mock()
        makeHeterogenousStub = mock()
        stub(RWRtoolkit::RWR_make_multiplex, 'make_homogenous_network', makeHomogenousStub)
        stub(RWRtoolkit::RWR_make_multiplex, 'make_heterogeneous_multiplex', makeHeterogenousStub)
        
        expect_error(RWRtoolkit::RWR_make_multiplex(flistFilePath))
        expect_called(makeHomogenousStub, 0)
        expect_called(makeHeterogenousStub, 0)
    })

    it('takes flist and makes a heterogeneous multiplex with default parameters', {
        ## Heterogeneous networks
        nw.groupsInput <- list_of( nwTibble_Het1,nwTibble_Het2, nwTibble_Het3 )
        flistFilePath <- '../testFlists/testFlist_heterogeneous.txt'

        makeHomogenousStub = mock()
        makeHeterogenousStub = mock()
        stub(RWRtoolkit::RWR_make_multiplex, 'make_homogenous_network', makeHomogenousStub)
        stub(RWRtoolkit::RWR_make_multiplex, 'make_heterogeneous_multiplex', makeHeterogenousStub)
        
        invisible(RWRtoolkit::RWR_make_multiplex(flistFilePath))

        expect_called(makeHomogenousStub, 0)
        expect_called(makeHeterogenousStub, 1)
        expect_args(makeHeterogenousStub, 1, nw.groupsInput, 0.5, 0.5, 'network.Rdata', FALSE)
    })
    
})


# removes test files so as not to clutter everything up
teardown({
    system('rm network.Rdata')
    system('rm testthatOutput.txt')
    system('rm testthatOutputP5.txt')
    system('rm testthatOutputP7.txt')
    system('rm testthatOutput1P0.txt')
    system('rm testthatOutputP5_fromTabDelmited.txt')
    system('rm testthatOutputP5_fromNoHeader.txt')
    system('rm testthatOutputP5_fromMixedDelim.txt')
    system('rm testthatOutputP5_fromMixedHeader.txt')
    system('rm testthatOutputP5_fromMixedHeaderDelim.txt')

}, env=parent.frame())

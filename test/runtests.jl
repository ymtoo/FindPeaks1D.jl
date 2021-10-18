using FindPeaks1D

using Test

@testset "FindPeaks1D" begin
    x = [13, 12, 14, 18, 19, 19, 19, 15, 11, 6, 4, 10, 8, 13, 8, 13, 3, 18, 7, 4]

    pkindices, _ = FindPeaks1D.localmaxima1d(x)
    @test pkindices == [6, 12, 14, 16, 18]

    pkindicesh1, _ = findpeaks1d(x; height=11)
    pkindicesh2, _ = findpeaks1d(x; height=(11, 16))
    @test pkindicesh1 == [6, 14, 16, 18]
    @test pkindicesh2 == [14, 16]

    pkindicesd1, _ = findpeaks1d(x, distance=1)
    pkindicesd2, _ = findpeaks1d(x, distance=3)
    @test pkindicesd1 == [6, 12, 14, 16, 18]
    @test pkindicesd2 == [6, 14, 18]

    prominences1, leftbases1, rightbases1 = peakprominences1d(x, pkindices)
    prominences2, leftbases2, rightbases2 = peakprominences1d(x, pkindices, 2)
    @test prominences1 == [7, 2, 9, 9, 14]
    @test leftbases1 == [2, 11, 11, 11, 17]
    @test rightbases1 == [17, 13, 17, 17, 20]
    @test prominences2 == [0, 2, 5, 5, 11]
    @test leftbases2 == [6, 11, 13, 15, 17]
    @test rightbases2 == [6, 13, 15, 17, 19]

    widths1, widthheight1, leftips1, rightips1 = peakwidths1d(x, pkindices, 1.0)
    @test all(isapprox.(widths1 ,[6.75, 1.3333, 5.900, 5.900, 2.9333], atol=0.001))
    @test widthheight1 == [12, 8, 4, 4, 4]
    @test all(isapprox.(leftips1 ,[2.0, 11.6667, 11.0, 11.0, 17.0667], atol=0.001))
    @test all(isapprox.(rightips1 ,[8.75, 13.0, 16.9, 16.9, 20.0], atol=0.001))

    pkindices1, _ = findpeaks1d(x, width=3.0, relheight=1.0)
    @test pkindices1 == [6,14,16]
    pkindices2, _ = findpeaks1d(x, width=2.0, relheight=0.5)
    @test pkindices2 == [6]
    pkindices3, _ = findpeaks1d(x, prominence=3)
    @test pkindices3 == [6,14,16,18]
    pkindices4, _ = findpeaks1d(x, prominence=100)
    @test pkindices4 == Int64[] 
    # test for less strict type signature
    pkindices1, _ = findpeaks1d(x, width=3//1, relheight=1.0)
    @test pkindices1 == [6,14,16]
    pkindices2, _ = findpeaks1d(x, width=2, relheight=1//2)
    @test pkindices2 == [6]
    pkindices3, _ = findpeaks1d(x, prominence=3//1)
    @test pkindices3 == [6,14,16,18]
    pkindices4, _ = findpeaks1d(x, prominence=100//1)
    @test pkindices4 == Int64[] 


    @test_throws ArgumentError pkindiceswlenerror, _ = findpeaks1d(x, prominence=0.1, wlen=0)

    x = [0, 1, 1, 1, 1, 1, 0]
    pkindices, _ = FindPeaks1D.localmaxima1d(x)
    @test length(pkindices) == 1
    @test pkindices[1] == 4

end

function tests = runTest
tests = functiontests(localfunctions);
close all;
end


function testRefrIndex(testCase)
act = refrindx(1064e-9,'air');
exp = 1.000274175520543;
verifyEqual(testCase,act,exp,'AbsTol',100*eps)

act = refrindx(1064e-9,'vacuum');
exp = 1;
verifyEqual(testCase,act,exp,'AbsTol',100*eps)

act = refrindx(1064e-9,'N-BAK 4');
exp = 1.556824316527775;
verifyEqual(testCase,act,exp,'AbsTol',100*eps)

end


function testExample1(testCase)
act = example1();
exp = 97.959183673469383;  % N-BAK 4
%exp = 1.028571428571429e+02;  % bk7
verifyEqual(testCase,act,exp,'AbsTol',100*eps(exp))

end

function testExample2(testCase)
actSolution = example2();
expSolution = 0.003573243838131;
verifyEqual(testCase,actSolution,expSolution,'AbsTol',expSolution/1e4)
end

function testExample2b(testCase)
actSolution = example2b();
expSolution = 1.137225471691192e-04;
verifyEqual(testCase,actSolution,expSolution,'AbsTol',expSolution/1e4)

end


function testExample3(testCase)
[act1,act2,act3] = example3();
exp1=[100,117.210229753348,137.382379588326,161.026202756094,188.739182213510,221.221629107045,259.294379740467,303.919538231320,356.224789026244,417.531893656040,489.390091847749,573.615251044868,672.335753649934,788.046281566991,923.670857187387,1082.63673387405,1268.96100316792,1487.35210729351,1743.32882219999,2043.35971785694,2395.02661998749,2807.21620394118,3290.34456231267,3856.62042116347,4520.35365636024,5298.31690628371,6210.16941891562,7278.95384398315,8531.67852417281,10000];
exp2=[9.34451105007655;9.49182955557099;9.63013637176348;9.75482897966147;9.86925000442919;9.97165940935692;10.0638066634026;10.1464060402606;10.2208006855949;10.2854728842908;10.3434459727582;10.3943171114229;10.4387954796068;10.4775595442909;10.5112487023554;10.5392430363891;10.5632189231423;10.5850055736504;10.6037843842852;10.6199473884088;10.6338416893220;10.6457729252669;10.6560089803338;10.6647837171887;10.6716397669617;10.6780734758335;10.6835789011381;10.6882879428817;10.6923144751726;10.7158708438953];
exp3=[0.00328604807095230;0.00323830311837706;0.00341805300030108;0.00337929626879949;0.00341208190849644;0.00335151714322298;0.00329588265687466;0.00325880677876516;0.00329110657300590;0.00323753365970787;0.00325653859686002;0.00327133861566085;0.00328292733195737;0.00329205386889041;0.00329928393526338;0.00324014579490953;0.00317405714755687;0.00317770677593686;0.00318066810785930;0.00318308351134040;0.00318506328107119;0.00318669335769910;0.00318804110706117;0.00318915965841837;0.00315572648526303;0.00315650138977597;0.00315715053630977;0.00315769564656743;0.00315815436777012;0.00316070402830657];

verifyEqual(testCase,act1,exp1,'AbsTol',100*eps(exp1))
verifyEqual(testCase,act2,exp2,'AbsTol',100*eps(exp2))
verifyEqual(testCase,act3,exp3,'AbsTol',100*eps(exp3))
end







function testExamples4(testCase)
% Compare screen images

act=uint16(example4());
%imwrite(uint16(act), 'tests/example04.png');
exp= imread('tests/example04.png');

verifyEqual(testCase,act,exp);

close all;
end

function testExamples5(testCase)
% No results, just check it works
example5();close all;
end

function testExamples6(testCase)
% No results, just check it works
example6();close all;
end

function testExamples7(testCase)
% No results, just check it works
example7();close all;
end

function testExamples8(testCase)
% No results, just check it works
example8();close all;
end
function testExamples9(testCase)
% No results, just check it works
example9();close all;
end

function testExamples10(testCase)
% No results, just check it works
example10();close all;
end

function testExamples11(testCase)
% No results, just check it works
example11();close all;
end

function testExamples12(testCase)
% No results, just check it works
example12();close all;
end

function testExamples13(testCase)
% No results, just check it works
example13();close all;
end

function testExamples14(testCase)
% No results, just check it works
example14();close all;
end

function testExamples15(testCase)
% No results, just check it works
example15();close all;
end

function testExamples16(testCase)
% No results, just check it works
example16();close all;
end

function testExamples17(testCase)
% No results, just check it works
example17();close all;
end

function testExamples18(testCase)
% No results, just check it works
example18();close all;
end

function testExamples19(testCase)
% No results, just check it works
example19();close all;
end
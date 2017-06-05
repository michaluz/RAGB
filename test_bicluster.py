from unittest import TestCase
from compute_bicliques import Bicluster
from compute_bicliques import Gene
import compute_bicliques
from bitarray import bitarray

__author__ = 'user'


class TestBicluster(TestCase):
  def test_contain(self):
    bic1 = Bicluster()
    bic1.addNode(1)
    bic1.addNode(2)
    bitarray1 = bitarray()
    bitarray1.append(1)
    bitarray1.append(1)
    bitarray1.append(0)
    bitarray1.append(1)
    bic1.setBitArray(bitarray1)

    bic2 = Bicluster()
    bic2.addNode(1)
    bic2.addNode(2)
    bic2.setBitArray(bitarray1)

    bic3 = Bicluster()
    bitarray2 = bitarray()
    bitarray2.append(1)
    bitarray2.append(1)
    bitarray2.append(1)
    bitarray2.append(1)
    bic3.addNode(1)
    bic3.addNode(2)
    bic3.addNode(3)
    bic3.setBitArray(bitarray2)

    bic4 = Bicluster()
    bic4.addNode(1)
    bic4.addNode(2)
    bic4.addNode(3)
    bic4.setBitArray(bitarray1)

    bic5 = Bicluster()
    bic5.addNode(1)
    bic5.addNode(2)
    bic5.setBitArray(bitarray2)

    #test for equals biclusters
    self.assertEqual(bic1.contain(bic2),1)
    #bic3 contains bic1 nodes and bitarray
    self.assertEqual(bic1.contain(bic3),1)
    #bic1 not contain bic3
    self.assertEqual(bic3.contain(bic1),0)
    #bic4 contain only nodes but not bitarray
    self.assertEqual(bic1.contain(bic4),1)
    #bic5 contaim only bitarray
    self.assertEqual(bic1.contain(bic5),1)

  def test_containInArray(self):
    bic1 = Bicluster()
    bic1.addNode(1)
    bic1.addNode(2)
    bitarray1 = bitarray()
    bitarray1.append(1)
    bitarray1.append(1)
    bitarray1.append(0)
    bitarray1.append(1)
    bic1.setBitArray(bitarray1)

    bic2 = Bicluster()
    bic2.addNode(1)
    bic2.addNode(2)
    bic2.setBitArray(bitarray1)

    bic3 = Bicluster()
    bitarray2 = bitarray()
    bitarray2.append(1)
    bitarray2.append(1)
    bitarray2.append(1)
    bitarray2.append(1)
    bic3.addNode(1)
    bic3.addNode(2)
    bic3.addNode(3)
    bic3.setBitArray(bitarray2)

    bic4 = Bicluster()
    bic4.addNode(1)
    bic4.addNode(2)
    bic4.addNode(3)
    bic4.setBitArray(bitarray1)

    bic5 = Bicluster()
    bic5.addNode(1)
    bic5.addNode(2)
    bic5.setBitArray(bitarray2)

    BicArray = []
    BicArray.append(bic1)
    BicArray.append(bic2)
    self.assertEqual(bic4.containInArray(BicArray),0)
    BicArray.append(bic3)
    self.assertEqual(bic4.containInArray(BicArray),1)

  def test_calculat_subsets(self):
      a = compute_bicliques.calculate_subsets(50)
      self.assertEqual((1,2,3,4) in a,True)
      self.assertEqual((17,21,22,23) in a,True)



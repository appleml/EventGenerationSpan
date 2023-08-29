# import math
# precision_list = [1,2,3,4,5]
# recall_list = [2,2, 3,4,5]
# F1_list = [5, 5,3,2,1]
#
# avg_prec = sum(precision_list)/len(precision_list)
# avg_recall = sum(recall_list)/len(recall_list)
# avg_f1 = sum(F1_list)/len(F1_list)
#
# print(avg_prec, "---", avg_recall, '---', avg_f1)
#
# square = list(map(lambda x: math.pow((x-avg_prec),2), precision_list))
# print(square)
# standard_deviation = math.sqrt(sum(square)/len(precision_list))
# print(standard_deviation)

listA = ["A", "B", "C", "A", "D", "B"]
listA.remove("B")
print(listA)

BMI <-  data.frame(
  gender = c("Male", "Female","Female"),
  height = c(152, 150, 165),
  weight = c(81,70, 78),
  Age = c(33,22,57)
)

print(BMI)


df <-   data.frame(
  letter1 = letters[1:20],
  letter2 = LETTERS[1:20],
  num1 = 1:20,
  num2 = 21:40,
  num3 = 41:60
)

View(df)
head(df, 10)
tail(df, 10)

class(BMI)
colnames(BMI)
rownames(BMI)
attributes(BMI)
dim(BMI)
nrow(BMI)
ncol(BMI)
length(BMI)

BMI

BMI[1,1]

BMI[2,1] <- "Male"
BMI[4,] <- list("Female",140, 72,68 ) # we use a list ... Why ?
BMI[,5] <- c(0,1,0,1)
colnames(BMI)
colnames(BMI) <- c("gender", "height", "weight", "Age", "smoker")
colnames(BMI)[4] <- "age"
rownames(BMI) <- c("Patient1","Patient2","Patient3","Patient4")


